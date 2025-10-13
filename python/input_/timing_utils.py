import time
import textwrap
from datetime import datetime

def start_timers():
    """
    Start CPU and wall-clock timers.

    Returns:
        tuple: (start_cpu_time, start_wall_time) — Initial CPU and wall-clock times.
    """    
    return time.process_time(), time.time()


def get_elapsed_times(start_cpu_time, start_wall_time, fortran_cpu_time):
    """
    Calculate total CPU and wall-clock times, including both Python and Fortran usage.

    Parameters:
        start_cpu_time (float): CPU time at the start of the Python process.
        start_wall_time (float): Wall-clock time at the start.
        fortran_cpu_time (float): CPU time used by the Fortran code.

    Returns:
        tuple: (cpu_h, cpu_m, cpu_s, wall_h, wall_m, wall_s) — Formatted CPU and wall times.
    """    
    # Calc final times
    end_cpu_time = (time.process_time() - start_cpu_time) + fortran_cpu_time
    end_wall_time = time.time() - start_wall_time

    # round and take h, m, sec
    cpu_h, remainder = divmod(round(end_cpu_time), 3600)
    cpu_m, cpu_s = divmod(remainder, 60)
    
    wall_h, remainder = divmod(round(end_wall_time), 3600)
    wall_m, wall_s = divmod(remainder, 60)

    return cpu_h, cpu_m, cpu_s, wall_h, wall_m, wall_s


def get_termination_message(success, stderr):
    """
    Determine the termination message based on the execution result.

    Parameters:
        success (bool): True if the code completed successfully, False otherwise.
        stderr (str): Standard error output from the execution.

    Returns:
        str: Termination message.
    """    
    if success:
        return "🎸 Ob-La-Di, Ob-La-Done! 🎶", "Normal Termination"
    else:
        if "Segmentation fault" in stderr or "memory allocation issue" in stderr:
            return "🚪 Knock, knock, knockin’ on debug’s door. 🛠️🛠️", "Error Termination"
        else:
            return "💀 Ob-La-Di, Ob-La-Doom! 💀", "Error Termination"

def print_execution_summary(code, start_cpu_time, start_wall_time, fortran_cpu_time, success, stderr, output_file, errors=None, citations=None):
    """
    Write an execution summary including timing information and termination message.

    Parameters:
        code (str): Name of the executed code.
        start_cpu_time (float): Python CPU time at start.
        start_wall_time (float): Wall-clock time at start.
        fortran_cpu_time (float): CPU time used by the Fortran code.
        success (bool): True if execution completed successfully, False otherwise.
        stderr (str): Standard error output from execution.
        output_file (str): Path to the output file.
        errors (list, optional): List of errors to include in the summary.
        citations (list, optional): List of citations to include in the summary.
    """    

    # Formatted times
    cpu_h, cpu_m, cpu_s, wall_h, wall_m, wall_s = get_elapsed_times(start_cpu_time, start_wall_time, fortran_cpu_time)

    termination_text, termination_status = get_termination_message(success, stderr)

    # Errors 
    sticks = " " + "-" * 80
    output_summary = []                                                               
    errors = errors or []                                                             
    if errors:
        output_summary.append("\n" + "Input Errors Detected".center(81) + "\n")
        output_summary.append(sticks)
        for error in errors:
            if "\n" in error:
                for line in error.splitlines():
                    wrapped_lines = textwrap.wrap(line, width=80, break_long_words=False)
                    for wl in wrapped_lines:
                        output_summary.append(" " + wl)
            else:
                wrapped_error = textwrap.wrap(error, width=80, break_long_words=False)  # Avvolge il testo a max 80 caratteri per riga
                for line in wrapped_error:
                    output_summary.append(" " + line)  # Aggiunge un singolo spazio a ogni riga
        output_summary.append(sticks)

    cpu_line  = f"    CPU Time: {str(cpu_h).rjust(5)} h {str(cpu_m).zfill(2)} min {str(cpu_s).zfill(2)} sec"
    wall_line = f"Elapsed Time: {str(wall_h).rjust(5)} h {str(wall_m).zfill(2)} min {str(wall_s).zfill(2)} sec"

    if citations and termination_status != "Error Termination":
        output_summary.append(" Required citations:\n")
        for c in citations:
            output_summary.append(c.rstrip('\r\n'))  
        output_summary.append(sticks)
    
    # Final output
    output_summary.extend([ 
        f"\n{termination_text.center(79)}\n", 
        sticks,
        cpu_line.rjust(81), 
        wall_line.rjust(81), 
        sticks, 
        f"{termination_status} of {code} in date {datetime.now().strftime('%d/%m/%Y at %H:%M:%S')}".rjust(81), 
        sticks ])
    
    # Output file
    with open(output_file, "a") as f:
        f.write("\n".join(output_summary) + "\n")
