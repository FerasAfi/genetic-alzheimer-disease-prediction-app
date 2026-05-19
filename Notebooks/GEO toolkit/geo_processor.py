import GEOparse
import pandas as pd
import numpy as np


class GEOProcessor:
    def __init__(self, file_path):
        self.file_path = file_path
        self.gse = None

    # --------------------
    # LOAD DATA
    # --------------------
    def load_data(self):
        self.gse = GEOparse.get_GEO(filepath=self.file_path)
        print(f"Loaded {self.gse.name}")
        return self.gse

    # --------------------
    # INFO
    # --------------------
    def show_info(self):
        print(
            f"The {self.gse.name} dataset contains {len(self.gse.gsms)} patients\n"
            f"Title: {self.gse.metadata['title'][0]}\n"
            f"Summary: {self.gse.metadata['summary'][0]}"
        )

    # --------------------
    # PATIENT IDS
    # --------------------
    def get_patient_ids(self):
        return list(self.gse.gsms.keys())

    # --------------------
    # PATIENT DATA
    # --------------------
    def show_patient_data(self, patient_id):
        print(f"Patient ID:\n{patient_id}\n")
        print(f"Patient gene data:\n{self.gse.gsms[patient_id].table}\n")
        print(f"Patient metadata:\n{self.gse.gsms[patient_id].metadata}")

    # --------------------
    # PATIENT METADATA
    # --------------------
    def show_patient_metadata(self, patient_id):
        metadata = self.gse.gsms[patient_id].metadata

        for key, value in metadata.items():
            print(f"{key}: {value}")

    # --------------------
    # COLUMNS
    # --------------------
    def get_table_columns(self):
        first_gsm = list(self.gse.gsms.values())[0]
        return first_gsm.table.columns.to_list()

    # --------------------
    # NORMALIZE
    # --------------------
    def normalize(self, col):
        return col.strip().lower().replace(" ", "_").replace("-", "_")

    # --------------------
    # EXPRESSION DATAFRAMES
    # --------------------
    def build_data_frames(self):

        columns = self.get_table_columns()
        index = columns[0]
        data_cols = columns[1:]  # skip index column
        dfs = {self.normalize(col): [] for col in data_cols}

        for gsm_name, gsm in self.gse.gsms.items():
            df = gsm.table.set_index(index)
            normalized_cols = { self.normalize(c): c for c in df.columns }

            for col in dfs:

                if col in normalized_cols:
                    real_col = normalized_cols[col]

                    series = df[real_col]
                    series.name = gsm_name

                    dfs[col].append(series)

        dfs = {
            col: pd.concat(series_list, axis=1)
            for col, series_list in dfs.items()
            if len(series_list) > 0
        }

        return dfs
    # --------------------
    # PATIENT LABELS
    # --------------------
    def get_patient_labels(self, patient_id):
        labels = []

        for key, value in self.gse.gsms[patient_id].metadata.items():
            if key.startswith("characteristics_ch"):
                labels.extend(value)

        print(f"\nPatient: {patient_id}")
        for label in labels:
            print(label)

        return patient_id, labels

    # --------------------
    # FULL METADATA TABLE
    # --------------------
    def build_gse_metadata(self):
        meta_rows = []

        for gsm_name, gsm in self.gse.gsms.items():
            values = {"patient_id": gsm_name}

            labels = []

            for key, val in gsm.metadata.items():
                if key.startswith("characteristics_ch"):
                    labels.extend(val)

            for label in labels:
                parts = label.split(":", 1)

                if len(parts) == 2:
                    key, value = parts
                    values[key.strip().lower()] = value.strip()
                else:
                    values[label.strip().lower()] = np.nan

            meta_rows.append(values)

        return pd.DataFrame(meta_rows)


    # -----------------------------
    # GET PLATFORM IDS (GPLs)
    # -----------------------------
    def get_platforms(self):
        return list(self.gse.gpls.keys())



    # -----------------------------
    # GET SPECIFIC GPL OBJECT
    # -----------------------------
    def get_platform(self, gpl_id):
        return self.gse.gpls[gpl_id]



    # -----------------------------
    # SHOW GPL TABLE COLUMNS
    # -----------------------------
    def get_platform_columns(self, gpl_id):
        gpl = self.get_platform(gpl_id)
        return gpl.table.columns.tolist()



    # -----------------------------
    # SHOW SELECTED COLUMNS
    # -----------------------------
    def show_platform_columns_data(self, gpl_id, columns, n=5):
        gpl = self.get_platform(gpl_id)
        return gpl.table[columns].head(n)



    # -----------------------------
    # BUILD PROBE -> GENE MAPPING
    # -----------------------------
    def build_probe_mapping(
        self,
        gpl_id,
        probe_col="ID",
        gene_col="Gene Symbol"
    ):

        gpl = self.get_platform(gpl_id)

        mapping = gpl.table[[probe_col, gene_col]].copy()

        probe_to_gene = dict(
            zip(mapping[probe_col], mapping[gene_col])
        )

        return probe_to_gene



    # -----------------------------
    # REPLACE PROBE IDS WITH GENES
    # -----------------------------
    def replace_probe_names(
        self,
        df,
        gpl_id,
        probe_col="ID",
        gene_col="Gene Symbol"
    ):

        probe_to_gene = self.build_probe_mapping(
            gpl_id,
            probe_col,
            gene_col
        )

        df = df.copy()

        df.index = df.index.map(probe_to_gene)


        return df