import subprocess


def run_gmx_cmd(arguments, prompt_input=None, print_output=True):
    try:
        result = subprocess.run(
            arguments,
            stdout=subprocess.PIPE,    # Capture stdout
            stderr=subprocess.STDOUT,  # Redirect stderr to stdout
            text=True,                 # Return strings instead of bytes
            input=prompt_input,
            check=True                 # Raise exception on non-zero exit
        )
        returncode, stdout = result.returncode, result.stdout
    except subprocess.CalledProcessError as e:
        # e.output contains the full error message (thanks to stderr=STDOUT)
        returncode, stdout = e.returncode, e.output

    if returncode != 0:
        raise RuntimeError(f'{" ".join(arguments[:2])} failed with return code {returncode}:\n{stdout}')
    else:
        if print_output:
            print(f"\n{stdout}")

    return returncode, stdout
