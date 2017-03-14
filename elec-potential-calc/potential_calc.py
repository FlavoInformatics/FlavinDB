import os
import subprocess

from contextlib import contextmanager

@contextmanager
def cd(new_dir):
    prev_dir = os.getcwd()
    os.chdir(os.path.expanduser(new_dir))
    try:
        yield
    finally:
        os.chdir(prev_dir)

def __get_dx_file(pdb_code, output_folder):
    """
    Run ABPS for a given protein and save the results to output_folder

    :param str pdb_code: four letter PDB protein identifier
    :param str output_folder: folder were we want to store resulting files.
        Files will be stored under <output_folder>/<pdb_code>/
    :return: None
    """

    print "--------------------------"
    print "Calculating potential for {}".format(pdb_code)

    results_folder_path = os.path.join(output_folder, pdb_code,)
    folder_path_FAD = os.path.join(results_folder_path, "FAD")
    folder_path_FMN = os.path.join(results_folder_path, "FMN")

    if not os.path.exists(results_folder_path):
        os.makedirs(results_folder_path)

    # TODO: create a PDB2PQRArgsBundle that will hold all pdb2pqr args info and will have a get_args() function to
    # populate the args list passed to the subprocess call, That way we don't have to hardcode args and it will be more
    # flexible for future use
    subprocess.call(["pdb2pqr", "--ff=amber", "--apbs-input", "--ph-calc-method=propka", "--ligand={}".format(results_folder_path), pdb_code,
                    os.path.join(results_folder_path, "{}.pqr".format(pdb_code))])


    


    with cd(results_folder_path):
        subprocess.call(["apbs", "{}.in".format(pdb_code)])


def get_potentials(pdb_codes, output_folder):
    for code in pdb_codes:
        __get_dx_file(code, output_folder)





