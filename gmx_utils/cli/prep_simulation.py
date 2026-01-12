import os
import sys
import glob
import time
import math
import shutil
import argparse
from general_utils import utils
from gmx_utils import gmx_utils
from gmx_utils import data

def initialize(args):
    parser = argparse.ArgumentParser(
        description="Prepare a GROMACS production simulation given necessary input files."
    )
    parser.add_argument(
        '-i',
        '--input_dir',
        type=str,
        required=True,
        help='The folder containing the topology (.top) file and initial structure (.gro) file.'
    )
    parser.add_argument(
        '-m',
        '--mdp_dir',
        type=str,
        help='The directory containing input MDP files for ion addition (ions.mdp), energy minimization (em.mdp), \
            NVT equilibration (nvt_equil.mdp), and NPT equilibration (npt_equil.mdp). If not specified, the default MDP files in \
            the package will be used.'
    )
    parser.add_argument(
        '-l',
        '--log',
        type=str,
        default='prep_simulation.log',
        help='Log file to record the output of the preparation steps. The default is prep_simulation.log.'
    )
    args_parse = parser.parse_args(args)

    return args_parse

def main():
    t1 = time.time()
    args = initialize(sys.argv[1:])
    sys.stdout = utils.Logger(args.log)
    sys.stderr = utils.Logger(args.log)

    print(f"\nCommand line: {' '.join(sys.argv)}")
    print(f"Current working directory: {os.getcwd()}")

    # Test GROMACS installation
    gmx_utils.run_gmx_cmd(['gmx', '--version'], print_output=False)

    # Step up directories and files
    mdp_dir = args.mdp_dir if args.mdp_dir else data.mdp_dir
    dir_list = ['topology', 'box', 'solv_ions', 'em', 'equil', 'production']
    for d in dir_list:
        os.makedirs(d, exist_ok=True)
    os.makedirs(os.path.join('equil', 'NVT'), exist_ok=True)
    os.makedirs(os.path.join('equil', 'NPT'), exist_ok=True)

    dst_topology = 'topology'
    for filename in os.listdir(args.input_dir):
        # Copy files from input_preparation/gmx_inputs to topology
        shutil.copy(os.path.join(args.input_dir, filename), os.path.join(dst_topology, filename))

    gro_list = glob.glob(os.path.join(dst_topology, '*.gro'))
    top_list = glob.glob(os.path.join(dst_topology, '*.top'))
    assert len(gro_list) == 1, "There should be exactly one GRO file in the input directory."
    assert len(top_list) == 1, "There should be exactly one TOP file in the input directory."

    input_gro = gro_list[0]
    input_top = top_list[0]

    print("\n1. Simulation box creation")
    print("==========================")
    gmx_args = [
        'gmx', 'editconf',
        '-f', input_gro,
        '-o', os.path.join('box', 'box.gro'),
        '-bt', 'cubic',
        '-d', '1.0',
        '-c'
    ]
    print(f'\nRunning command: {" ".join(gmx_args)}')
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args)
    box_volume = float(stdout.split(':')[-1].split()[0])

    print("2. Solvation")
    print("============")
    gmx_args = [
        'gmx', 'solvate',
        '-cp', os.path.join('box', 'box.gro'),
        '-o', os.path.join('solv_ions', 'solv.gro'),
        '-p', input_top,
        '-cs'
    ]
    print(f'\nRunning command: {" ".join(gmx_args)}')
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args)

    print("3. Neutralization with ions")
    print("============================")
    gmx_args = [
        'gmx', 'grompp',
        '-f', os.path.join(mdp_dir, 'ions.mdp'),
        '-c', os.path.join('solv_ions', 'solv.gro'),
        '-p', input_top,
        '-o', os.path.join('solv_ions', 'ions.tpr'),
        '-maxwarn', '1'
    ]
    print(f'\nRunning command 1: {" ".join(gmx_args)}')
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args)

    for line in stdout.splitlines():
        if "System has non-zero total charge:" in line:
            total_charge = float(line.split(":")[-1])
            break

    n_ions = int(math.ceil(box_volume * 1E-27 * 1000 * 0.15 * 6.022E23))
    if total_charge < 0:
        n_chloride = n_ions
        n_sodium = n_ions + int(abs(total_charge))
    else:
        n_sodium = n_ions
        n_chloride = n_ions + int(abs(total_charge))

    gmx_args = [
        'gmx', 'genion',
        '-s', os.path.join('solv_ions', 'ions.tpr'),
        '-o', os.path.join('solv_ions', 'ions.gro'),
        '-p', input_top,
        '-pname', 'NA',
        '-nname', 'CL',
        '-np', str(n_sodium),
        '-nn', str(n_chloride),
    ]
    print(f'\nRunning command 2: {" ".join(gmx_args)}')
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args, prompt_input='SOL\n')

    print("4. Energy minimization")
    print("======================")
    gmx_args = [
        'gmx', 'grompp',
        '-f', os.path.join(mdp_dir, 'em.mdp'),
        '-c', os.path.join('solv_ions', 'ions.gro'),
        '-p', input_top,
        '-o', os.path.join('em', 'em.tpr'),
        '-maxwarn', '1'
    ]
    print(f'\nRunning command 1: {" ".join(gmx_args)}')
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args)

    gmx_args = ['gmx', 'mdrun', '-deffnm', 'em/em']
    print(f'\nRunning command 2: {" ".join(gmx_args)}')
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args, prompt_input=None)

    print("5. NVT Equilibration")
    print("===================")
    gmx_args = [
        'gmx', 'grompp',
        '-f', os.path.join(mdp_dir, 'nvt_equil.mdp'),
        '-c', os.path.join('em', 'em.gro'),
        '-p', input_top,
        '-o', os.path.join('equil', 'NVT', 'equil.tpr'),
        '-maxwarn', '1'
    ]
    print(f'\nRunning command 1: {" ".join(gmx_args)}')
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args)

    gmx_args = ['gmx', 'mdrun', '-deffnm', 'equil/NVT/equil']
    print(f'\nRunning command 2: {" ".join(gmx_args)}')
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args, prompt_input=None)

    print("6. NPT Equilibration")
    print("===================")
    gmx_args = [
        'gmx', 'grompp',
        '-f', os.path.join(mdp_dir, 'npt_equil.mdp'),
        '-c', os.path.join('equil', 'NVT', 'equil.gro'),
        '-p', input_top,
        '-o', os.path.join('equil', 'NPT', 'equil.tpr'),
        '-maxwarn', '2'
    ]
    print(f'\nRunning command 1: {" ".join(gmx_args)}')
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args)

    gmx_args = ['gmx', 'mdrun', '-deffnm', 'equil/NPT/equil']
    print(f'\nRunning command 2: {" ".join(gmx_args)}')
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args, prompt_input=None)

    print(f"Elapsed time: {utils.format_time(time.time() - t1)}")
    print("Preparation completed successfully.")
