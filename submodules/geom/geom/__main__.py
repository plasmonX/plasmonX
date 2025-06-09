import sys 

from .classes import input_class
from .functions import general, translate, rotate, create_geom, various, output


#                    ██████╗ ███████╗ ██████╗ ███╗   ███╗     ██████╗ ██████╗ ██████╗ ███████╗                
#                   ██╔════╝ ██╔════╝██╔═══██╗████╗ ████║    ██╔════╝██╔═══██╗██╔══██╗██╔════╝                
#   █████╗█████╗    ██║  ███╗█████╗  ██║   ██║██╔████╔██║    ██║     ██║   ██║██║  ██║█████╗      █████╗█████╗
#   ╚════╝╚════╝    ██║   ██║██╔══╝  ██║   ██║██║╚██╔╝██║    ██║     ██║   ██║██║  ██║██╔══╝      ╚════╝╚════╝
#                   ╚██████╔╝███████╗╚██████╔╝██║ ╚═╝ ██║    ╚██████╗╚██████╔╝██████╔╝███████╗                
#                    ╚═════╝ ╚══════╝ ╚═════╝ ╚═╝     ╚═╝     ╚═════╝ ╚═════╝ ╚═════╝ ╚══════╝                
   

# ============================================================================================================ #
#                                        Program by Pablo Grobas Illobre                                       #
#                                                                                                              #
#                                       Contact: pgrobasillobre@gmail.com                                      #
# ============================================================================================================ #

# ---------------------------------------------------------
# PURPOSE: Create XYZ files of metal nanoparticles/graphene
#          structures and manage geometries
#
# EXECUTION details: python3 geom -h 
# ---------------------------------------------------------

def main():
    """
    Main function to initialize input parameters and execute the appropriate geometry processing task.

    Returns:
        None: Calls the relevant function based on the user's input.
    """
    try:
        # Initialize input class and parse command-line arguments
        inp = input_class.input_class()
        general.read_command_line(sys.argv, inp)

        # Select and execute the appropriate task
        if inp.translate:
            translate.select_case(inp)
        elif inp.rotate:
            rotate.select_case(inp)
        elif inp.create_geom:
            create_geom.select_case(inp)
        elif inp.small_tasks:
            various.select_case(inp)
        else:
            output.error("No valid task specified. Use -h for help.")

    except Exception as e:
        output.error(f"An error occurred: {e}")


if __name__ == "__main__":
    main()
