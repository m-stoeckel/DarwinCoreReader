import argparse
import csv
import os
import subprocess
from pathlib import Path
from typing import Dict, IO, List

from regex import regex
from tqdm import tqdm

kingdom_class_map = {
    "Animalia": "Animal_Fauna",
    "Archaea": "Archaea",
    "Bacteria": "Bacteria",
    "Chromista": "Chromista",
    "Fungi": "Fungi",
    "Lichen": "Lichen",
    "Plantae": "Plant_Flora",
    "Protozoa": "Protozoa",
    "Taxon": "Taxon",
    "Viruses": "Viruses"
}

class_kingdom_map = {
    "Animal_Fauna": "Animalia",
    "Archaea": "Archaea",
    "Bacteria": "Bacteria",
    "Chromista": "Chromista",
    "Fungi": "Fungi",
    "Lichen": None,
    "Plant_Flora": "Plantae",
    "Protozoa": "Protozoa",
    "Taxon": None,
    "Viruses": "Viruses"
}

vernacular_name_sepearator: regex.Pattern = regex.compile(", ?")
strip_author_pattern: regex.Pattern = regex.compile(r"^([\p{upper}]?[\P{upper}]+) (in)?.*$", flags=regex.UNICODE)


def linecount(filename):
    with open(filename, 'r') as fin:
        count = sum(1 for _ in fin)
    return count


class DarwinCoreProcessor:
    def __init__(self, filtered_taxon_tsv, vernacular_name_tsv, base_path_name, base_file_name,
                 separate_vernacular_names=True, use_subfolders=True, file_extension=".list", filter_languages=None,
                 sort=True, delete_empty=True):
        self.taxon_tsv = filtered_taxon_tsv
        self.vernacular_name_tsv = vernacular_name_tsv
        self.base_path: Path = Path(base_path_name)
        self.base_file_name = base_file_name
        self.use_subfolders = use_subfolders
        self.separate_vernacular_names = separate_vernacular_names
        self.file_extension = file_extension

        self.id_kingdom_map: Dict[str, str] = dict()
        if filter_languages is None:
            filter_languages = {"de", "la"}
        self.filter_languages = filter_languages

        if not self.separate_vernacular_names:
            self.out_files_tx = self.get_out_files("")
            self.out_files_vn = self.out_files_tx
        else:
            self.out_files_tx = self.get_out_files("_tx")
            self.out_files_vn = self.get_out_files("_vn")

        # fields to be defined in subclasses
        self.taxon_id = None
        self.canonical_name_id = None
        self.kingdom_id = None

        self.process()
        if sort:
            self.sort()
        if delete_empty:
            self.delete_empty()

    def get_out_files(self, identifier) -> Dict[str, IO]:
        out_files_dict: Dict[str, IO] = {}
        for clazz in class_kingdom_map.keys():
            if self.use_subfolders:
                clazz_folder = self.base_path / clazz
                clazz_folder.mkdir(exist_ok=True)
                file_path = clazz_folder / f"{self.base_file_name}{identifier}{self.file_extension}"
            else:
                file_path = self.base_path / f"{self.base_file_name}_{clazz.lower()}{identifier}{self.file_extension}"
            out_files_dict[clazz] = file_path.open("w", encoding="utf8")
        return out_files_dict

    def process(self):
        # taxa file processing
        lc = linecount(self.taxon_tsv)
        with open(self.taxon_tsv, 'r', encoding='utf8') as fin:
            reader = csv.reader(tqdm(fin, desc="Processing taxa file", total=lc), delimiter='\t',
                                quoting=csv.QUOTE_NONE)
            self.read_taxon_tsv(reader)

        # vernacular name file processing
        lc = linecount(self.vernacular_name_tsv)
        with open(self.vernacular_name_tsv, 'r', encoding='utf8') as fin:
            reader = csv.reader(tqdm(fin, desc="Processing vernacular name file", total=lc),
                                delimiter='\t', quoting=csv.QUOTE_NONE)
            self.read_vernacular_names_tsv(reader)

        for fout in self.out_files_tx.values():
            fout.close()
        if self.separate_vernacular_names:
            for fout in self.out_files_vn.values():
                fout.close()

    def read_taxon_tsv(self, reader: csv.reader):
        raise NotImplementedError()

    def process_taxon_header(self, reader):
        raise NotImplementedError()

    def get_taxon_id(self, row):
        return row[self.taxon_id]

    def get_name(self, row):
        return row[self.canonical_name_id]

    def get_kingdom(self, row: List[str], name: str):
        """
        Get the kingdom of this taxon. Categorize all taxa containing the word "Lichen" as a member of this family.

        :param row: The current row.
        :param name: The taxon's name
        :return:
        """
        return row[self.kingdom_id]

    def read_vernacular_names_tsv(self, reader: csv.reader):
        raise NotImplementedError

    def process_vernacular_names_header(self, reader: csv.reader):
        for header in reader:
            self.taxon_id = int(header.index('taxonID'))
            self.language_id = int(header.index('language'))
            self.vernacular_name_id = int(header.index('vernacularName'))
            break

    def get_vernacular_name(self, row):
        return row[self.vernacular_name_id]

    def get_language(self, row):
        return row[self.language_id]

    def sort(self):
        file_names = [file.name for file in self.out_files_tx.values()]
        if self.separate_vernacular_names:
            file_names += [file.name for file in self.out_files_vn.values()]
        for file_name in tqdm(file_names, desc="Sorting files"):
            subprocess.call(f"sort -u -o {file_name} {file_name}", shell=True)

    def delete_empty(self):
        file_names = [file.name for file in self.out_files_tx.values()]
        if self.separate_vernacular_names:
            file_names += [file.name for file in self.out_files_vn.values()]
        for file_name in tqdm(file_names, desc="Deleting empty files"):
            if os.stat(file_name).st_size == 0:
                os.remove(file_name)


class BackboneProcessor(DarwinCoreProcessor):
    base_uri = "https://www.gbif.org/species/"

    # table header keys for the GBIF Backbone Taxonomy
    taxon_id_key = 'taxonID'
    kingdom_key = 'kingdom'
    canonical_name_key = 'canonicalName'
    generic_name_key = 'genericName'

    def __init__(self, filtered_taxon_tsv, vernacular_name_tsv, base_path_name, base_file_name="gbif_backbone",
                 **kwargs):
        super().__init__(filtered_taxon_tsv, vernacular_name_tsv, base_path_name, base_file_name, **kwargs)

    def read_taxon_tsv(self, reader: csv.reader):
        self.process_taxon_header(reader)

        for line in reader:
            taxon_id = self.get_taxon_id(line)
            name = self.get_name(line)
            if name is None or name == "" or name.isspace():
                continue

            # get the kingdom with a Lichen workaround
            kingdom = self.get_kingdom(line, name)

            # save the taxon id for vernacular name mapping
            self.id_kingdom_map[taxon_id] = kingdom

            # get the TTLab class for this kingdom, defaulting to Taxon
            clazz = kingdom_class_map.get(kingdom, "Taxon")

            self.out_files_tx[clazz].write(f"{name}\t{self.base_uri}{taxon_id}\n")

    def process_taxon_header(self, reader):
        for header in reader:
            self.taxon_id = int(header.index(self.taxon_id_key))
            self.canonical_name_id = int(header.index(self.canonical_name_key))
            self.kingdom_id = int(header.index(self.kingdom_key))
            break

    def process_vernacular_names_header(self, reader: csv.reader):
        for header in reader:
            self.taxon_id = int(header.index('taxonID'))
            self.language_id = int(header.index('language'))
            self.vernacular_name_id = int(header.index('vernacularName'))
            break

    def read_vernacular_names_tsv(self, reader: csv.reader):
        self.process_vernacular_names_header(reader)

        filtered = 0
        all_count = 0
        for row in reader:
            taxon_id = self.get_taxon_id(row)
            vernacular_name = self.get_vernacular_name(row)
            lang = self.get_language(row)
            if any((lang == filter_lang) for filter_lang in self.filter_languages):
                for name in vernacular_name_sepearator.split(vernacular_name):
                    clazz = kingdom_class_map.get(self.id_kingdom_map.get(taxon_id, ""), "Taxon")
                    self.out_files_vn[clazz].write(f"{name}\t{self.base_uri}{taxon_id}\n")
            else:
                filtered += 1
            all_count += 1

        print(f"Filtered {filtered}/{all_count} vernacular names on selected languages: {self.filter_languages}")


class CatalogueOfLifeProcessor(DarwinCoreProcessor):
    base_uri = "http://www.catalogueoflife.org/col/details/species/id/"
    taxon_id_key = 'taxonID'
    identifier_id_key = 'identifier'
    accepted_name_usage_id_key = 'acceptedNameUsageID'
    kingdom_key = 'kingdom'
    canonical_name_key = 'canonicalName'
    scientific_name_key = 'scientificName'
    generic_name_key = 'genericName'
    genus_key = 'genus'
    specific_epithet_key = 'specificEpithet'
    infraspecific_epithet_key = 'infraspecificEpithet'

    def __init__(self, filtered_taxon_tsv, vernacular_name_tsv, base_path_name, base_file_name="catalogue_of_life",
                 build_from_generic_name=True, strip_author=False, filter_languages=None, **kwargs):
        self.build_from_generic_name = build_from_generic_name
        self.strip_author = strip_author
        self.taxon_id_identifier_mapping = {}
        self.accepted_name_taxon_mapping = {}
        if filter_languages is None:
            filter_languages = {"German", "Ger", "Prussian, Old DE", "Prussian"}
        super().__init__(filtered_taxon_tsv, vernacular_name_tsv, base_path_name, base_file_name,
                         filter_languages=filter_languages, **kwargs)

    def read_taxon_tsv(self, reader: csv.reader) -> None:
        self.process_taxon_header(reader)

        for row in reader:
            taxon_id = self.get_taxon_id(row)
            identifier = self.get_identifier(row)
            accepted_name_usage = row[self.accepted_name_usage_id]
            self.taxon_id_identifier_mapping[taxon_id] = identifier
            name = self.get_name(row)

            if name is None or name == "" or name.isspace() or name == "incertae sedis":
                continue

            # get the kingdom with a Lichen workaround
            kingdom = self.get_kingdom(row, name)

            # save the taxon id for vernacular name mapping
            self.id_kingdom_map[taxon_id] = kingdom

            # get the TTLab class for this kingdom, defaulting to Taxon
            clazz = kingdom_class_map.get(kingdom, "Taxon")

            if identifier == "" or identifier.isspace():
                other_taxa = self.accepted_name_taxon_mapping.get(accepted_name_usage, [])
                self.accepted_name_taxon_mapping.update({accepted_name_usage: other_taxa + [(name, taxon_id)]})
            else:
                self.out_files_tx[clazz].write(f"{name}\t{self.base_uri}{identifier}\n")
                for other_name, other_taxon_id in self.accepted_name_taxon_mapping.get(taxon_id, []):
                    self.out_files_tx[clazz].write(f"{other_name}\t{self.base_uri}{identifier}\n")
                    self.taxon_id_identifier_mapping[other_taxon_id] = identifier

    def process_taxon_header(self, reader) -> None:
        for header in reader:
            self.taxon_id = int(header.index(self.taxon_id_key))
            self.identifier_id = int(header.index(self.identifier_id_key))
            self.accepted_name_usage_id = int(header.index(self.accepted_name_usage_id_key))
            self.scientific_name_id = int(header.index(self.scientific_name_key))
            self.kingdom_id = int(header.index(self.kingdom_key))

            self.generic_name_id = int(header.index(self.generic_name_key))
            self.genus_id = int(header.index(self.genus_key))  # TODO: Add a genus + sepi + isepi?
            self.specific_epithet_id = int(header.index(self.specific_epithet_key))
            self.infraspecific_epithet_id = int(header.index(self.infraspecific_epithet_key))
            break

    def get_name(self, row) -> str:
        if self.build_from_generic_name:
            generic_name = row[self.generic_name_id]
            specific_epithet = row[self.specific_epithet_id]
            infraspecific_epithet = row[self.infraspecific_epithet_id]
            name = " ".join([generic_name, specific_epithet, infraspecific_epithet])
        else:
            name = row[self.scientific_name_id]
            if self.strip_author:
                match = strip_author_pattern.match(name)
                if match is not None:
                    name = match.group(1)
        return name

    def get_identifier(self, row):
        return row[self.identifier_id]

    def read_vernacular_names_tsv(self, reader: csv.reader):
        self.process_vernacular_names_header(reader)

        filtered = 0
        all_count = 0
        for row in reader:
            taxon_id = self.get_taxon_id(row)
            identifier = self.taxon_id_identifier_mapping.get(taxon_id, taxon_id)
            vernacular_name = self.get_vernacular_name(row)
            lang = self.get_language(row)
            if any((lang == filter_lang) for filter_lang in self.filter_languages):
                for name in vernacular_name_sepearator.split(vernacular_name):
                    clazz = kingdom_class_map.get(self.id_kingdom_map.get(taxon_id, ""), "Taxon")
                    # if identifier == "" or identifier.isspace:
                    #     self.out_files_vn[clazz].write(f"{name}\t{self.gbif_uri}{taxon_id}\n")
                    # else:
                    self.out_files_vn[clazz].write(f"{name}\t{self.base_uri}{identifier}\n")
            else:
                filtered += 1
            all_count += 1

        print(f"Filtered {filtered}/{all_count} vernacular names on selected languages: {self.filter_languages}")


if __name__ == '__main__':
    BackboneProcessor("/home/manu/Work/data/datasets/backbone-current/filtered-Taxon.tsv",
                      "/home/manu/Work/data/datasets/backbone-current/filtered-VernacularName.tsv",
                      "/home/manu/Work/data/classes/biofid/")
    CatalogueOfLifeProcessor("/home/manu/Work/data/datasets/archive-complete/taxa.txt",
                             "/home/manu/Work/data/datasets/archive-complete/vernacular.txt",
                             "/home/manu/Work/data/classes/biofid/")
    exit()
    parser = argparse.ArgumentParser()
