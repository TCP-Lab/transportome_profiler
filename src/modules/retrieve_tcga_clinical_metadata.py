"""Ask the GDC data portal for clinical information. Unpack it (clean it) from
the JSON and paste it inside a Pandas Dataframe together with the standard
patient identifier.
"""

import json
import logging
import time
from functools import reduce

import pandas as pd
import requests
from sys import stdout

logging.basicConfig(format="%(asctime)s %(funcName)s@%(filename)s [%(levelname)s]: %(message)s", level=logging.INFO)
log = logging.getLogger("asclepius") # A famous greek clinician

ALL_STUDIES = [f"TCGA-{x}" for x in (
        "LAML",
        "ACC",
        "BLCA",
        "LGG",
        "BRCA",
        "CESC",
        "CHOL",
        #"LCML",
        "COAD",
        "ESCA",
        "GBM",
        "HNSC",
        "KICH",
        "KIRC",
        "KIRP",
        "LIHC",
        "LUAD",
        "LUSC",
        "DLBC",
        "MESO",
        "OV", 
        "PAAD",
        "PCPG",
        "PRAD",
        "READ",
        "SARC",
        "SKCM",
        "STAD",
        "TGCT",
        "THYM",
        "THCA",
        "UCS",
        "UCEC",
        "UVM"
    )]

def call_portal(project_id: str, number: int = 1_000_000):
    """Retrieves clinical data regarding a project from the GDC database

    Requires an internet connection. Downloads data regarding all patients in
    the project (up to 'number'), and will then merge it into a single pandas
    DataFrame for manipulation.

    Writes all information in a log file in the same folder as the script.

    Args:
        project_id : The project's ID, such as TCGA-BRCA
        number: Download the first "number" of patients, up to all patients.
            This should generally be kept as a large enough values. Defaults
            to 1_000_000

    Returns:
        Pandas dataframe containing clinical information
    """

    cases_endpt = "https://api.gdc.cancer.gov/cases/"
    filters = {
        "op": "=",
        "content": {"field": "project.project_id", "value": project_id},
    }
    data_types = ["diagnoses", "demographic", "exposures"]

    dataframes = []

    for data_type in data_types:
        params = {
            "filters": json.dumps(filters),
            "format": "JSON",
            "expand": data_type,
            # Omission of the int makes the value have a dot (100.0)
            # and it causes an **internal server error**!!! 
            # WTF GDC, come on.
            "size": str(int(number)),
        }
        log.info(f"Calling GDC portal for {data_type} data")
        start = time.perf_counter()
        response = requests.get(cases_endpt, params=params)
        response.raise_for_status()
        log.info(
            f"Received {len(response.content)} bytes in {(time.perf_counter() - start):.2f} seconds"
        )
        warnings = json.loads(response.content.decode("utf-8"))["warnings"]
        if warnings:
            log.warning(f"There were some warnings when downloading data: {warnings}")
        # Decode and unpack the JSON response
        decoded_response = json.loads(response.content.decode("utf-8"))["data"]["hits"]
        # Clean up the data, and put it in a dataframe ------------------------
        missing_diagnoses = []
        cases = []
        for patient in decoded_response:
            clean_data = {}
            try:
                if data_type != "demographic":
                    # Diagnoses and exposures are in a list of 1 element, so
                    # I'm unlisting them here (the [0])
                    clean_data.update(patient[data_type][0])
                else:
                    # Demographic is just a dictionary, no need to unlist
                    clean_data.update(patient[data_type])
            except KeyError:
                missing_diagnoses.append(patient["submitter_id"])
            # Add the relevant patient ID to the cleaned data for merging
            clean_data.update({"submitter_id": patient["submitter_id"]})
            cases.append(clean_data)
        # Warn the user if something went wrong when retrieving the data
        if missing_diagnoses:
            str_missing_diagnoses = ", ".join(missing_diagnoses)
            log.warning(
                f"I found one or more missing {data_type}: {str_missing_diagnoses}"
            )
        # Finally, add the dataframe to the dataframe list
        dataframes.append(pd.DataFrame(cases))
    log.info("Dropping useless columns")
    clean_dframes = []
    for i in dataframes:
        clean_dframes.append(i.drop(columns = ["created_datetime", "state", "updated_datetime"]))
    dataframes = clean_dframes
    # Collapse all dataframes into a single one
    log.info("Collapsing received data")
    merged_frame = reduce(
        lambda x, y: pd.merge(x, y, on="submitter_id", how="outer"), dataframes
    )
    # Add the study ID
    merged_frame["study"] = project_id
    return merged_frame


def get_clinical_data(project_id: str, number: int):
    """Retrieves clinical data from the GDC data portal given a TCGA ID

    Gets data from the first NUMBER patients (default to 1 Million) in the
    TCGA project with id PROJECT_ID, and save them in csv format to OUTPUT_FILE.

    The missing values are replaced with the string "not reported",
    like TCGA does with their missing variables.
    """
    
    if project_id != "ALL":
        call_portal(project_id, number).to_csv(stdout, index = False)
        return
    
    all_data = []
    for id in ALL_STUDIES:
        log.info(f"Retrieving {id}..")
        all_data.append(call_portal(id, number))

    log.info("Collating all studies together")
    # use pd.concat
    
    merged_frame = pd.concat(all_data)
    

    merged_frame.to_csv(stdout, index=False)



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()

    parser.add_argument("project_id", help="TCGA shorthand project to download, e.g. 'TCGA-BRCA'. If 'ALL', downloads all studies.")
    parser.add_argument("--number", help="Number of patients to retrieve", default=1e6)

    args = parser.parse_args()

    get_clinical_data(args.project_id, args.number)

