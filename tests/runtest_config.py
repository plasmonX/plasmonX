def configure(options, input_files,extra_arg):
    """
    This function is used by runtest to configure runtest
    at runtime for code specific launch command and file naming.
    """

    from os import path
    from sys import platform

    launcher = 'plasmonX.py'
    launcher_full_path = path.normpath(path.join(options.binary_dir, launcher))

    if len(input_files) == 1:
        (inp) = input_files
        extra = None

    inp_no_prefix = inp[0]
    command = []
    command.append(launcher_full_path)
    command.append("-i "+inp_no_prefix)

    full_command = ' '.join(command)
    print(full_command)

    name_input_file = inp[0].split()[0]

    output_prefix = name_input_file[:-5]
    print(output_prefix)

    relative_reference_path = 'reference'

    return launcher, full_command, output_prefix, relative_reference_path

def configure_analysis(options, input_files,extra_arg):
    """
    This function is used by runtest to configure runtest
    at runtime for code specific launch command and file naming.
    """

    from os import path
    from sys import platform

    launcher = 'plasmonX_analysis.py'
    launcher_full_path = path.normpath(path.join(options.binary_dir, launcher))

    if len(input_files) == 1:
        (inp) = input_files
        extra = None

    inp_no_prefix = inp[0]
    command = []
    command.append(launcher_full_path)
    command.append("-i "+inp_no_prefix)

    full_command = ' '.join(command)
    print(full_command)

    name_input_file = inp[0].split()[0]
    output_prefix = "post_process_"+name_input_file[:-7]+"/"#+name_input_file[:-7]
    print(output_prefix)

    relative_reference_path = 'reference'

    return launcher, full_command, output_prefix, relative_reference_path
