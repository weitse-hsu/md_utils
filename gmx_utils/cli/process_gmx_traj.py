import os
import sys
import time
import argparse
from general_utils import utils
from gmx_utils import gmx_utils


def initialize(args):
    parser = argparse.ArgumentParser(
        description="Process a GROMACS trajectory by recentering and rewrapping molecules."
    )
    parser.add_argument(
        '-i',
        '--input',
        type=str,
        required=True,
        help='Input GROMACS trajectory file.'
    )
    parser.add_argument(
        '-t',
        '--tpr',
        type=str,
        help='Input GROMACS TPR file corresponding to the trajectory. If not provided, it is assumed to \
            have the same prefix as the trajectory file with .tpr extension.'
    )
    parser.add_argument(
        '-o',
        '--output',
        type=str,
        help='Output GROMACS trajectory file. The default is the original filename with _center.pdb suffix.'
    )
    parser.add_argument(
        '-g',
        '--grps',
        type=str,
        nargs=4,
        help='Four index groups for centering and output selection for the two trjconv commands. The default \
            is "Backbone System Backbone System".'
    )
    parser.add_argument(
        '-dt',
        '--time_step',
        type=float,
        default=200,
        help='Time step between frames in the trajectory (in ps). If not provided, the time step is read \
            from the trajectory file. The default is 200 ps.'
    )
    parser.add_argument(
        '-l',
        '--log',
        type=str,
        default='process_gmx_traj.log',
        help='Log file to record the output of the processing steps. The default is process_gmx_traj.log.'
    )
    args = parser.parse_args(args)
    return args


def main():
    t1 = time.time()
    args = initialize(sys.argv[1:])
    sys.stdout = utils.Logger(args.log)
    sys.stderr = utils.Logger(args.log)

    print(f"\nCommand line: {' '.join(sys.argv)}")
    print(f"Current working directory: {os.getcwd()}")

    # Test GROMACS installation
    gmx_utils.run_gmx_cmd(['gmx', '--version'], print_output=False)

    input_traj = args.input
    prefix = input_traj.rsplit('.', 1)[0]
    tpr_file = args.tpr if args.tpr else f"{prefix}.tpr"
    output_traj = args.output if args.output else f"{prefix}_center.xtc"
    grps = args.grps if args.grps else ['Backbone', 'System', 'Backbone', 'System']
    cmd_list = []

    # Step 1: Recenter and remove jumps
    gmx_args = [
        'gmx', 'trjconv',
        '-s', tpr_file,
        '-f', input_traj,
        '-o', f'{prefix}_nojump.xtc',
        '-center',
        '-pbc', 'nojump'
    ]

    if args.time_step:
        gmx_args.extend(['-dt', str(args.time_step)])

    print(f'\nRunning command: {" ".join(gmx_args)}')
    cmd_list.append(' '.join(gmx_args))
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args, prompt_input=f'{grps[0]}\n{grps[1]}\n')

    gmx_args = [
        'gmx', 'trjconv',
        '-s', tpr_file,
        '-f', f'{prefix}_nojump.xtc',
        '-o', output_traj,
        '-center',
        '-pbc', 'whole',
        '-ur', 'compact',
    ]

    if args.time_step:
        gmx_args.extend(['-dt', str(args.time_step)])

    print(f'\nRunning command: {" ".join(gmx_args)}')
    cmd_list.append(' '.join(gmx_args))
    returncode, stdout = gmx_utils.run_gmx_cmd(gmx_args, prompt_input=f'{grps[2]}\n{grps[3]}\n')

    print("Summary of commands executed:")
    for cmd in cmd_list:
        print(cmd)

    print(f"Elapsed time: {utils.format_time(time.time() - t1)}")
    print("Trajectory processing completed successfully.")
