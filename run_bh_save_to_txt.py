import subprocess
import shutil
import time
import os
import re

def remove_first_three_lines(file_path):
    # Read the file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Skip the first three lines and write the rest to the file
    with open(file_path, 'w') as file:
        file.writelines(lines[3:]) 

def move_text_file(source_folder, destination_folder, file_name):
    # Create the full source and destination paths
    source_path = os.path.join(source_folder, file_name)
    destination_path = os.path.join(destination_folder, file_name)

    # Check if the source file exists
    if not os.path.exists(source_path):
        print(f"Error: {file_name} not found in the source folder.")
        return
    
    # Check if the destination folder exists, if not, create it
    if not os.path.exists(destination_folder):
        os.makedirs(destination_folder)
        print(f"Destination folder {destination_folder} created.")
    
    # Move the file
    try:
        shutil.move(source_path, destination_path)
        print(f"Moved {file_name} to {destination_folder}")
    except Exception as e:
        print(f"Error moving {file_name}: {e}")


def execute_BH(working_directory, parameter_file, mass, dof):
    try:
        # Run the command with parameters.txt as an argument
        result = subprocess.run(["./BlackHawk_tot.x", "parameters.txt"],cwd = current_dir, text=True, capture_output=True, check=True)
        
        # Print the output and errors (if any)
        print("Output:\n", result.stdout)
        print("Errors:\n", result.stderr)
    
    except subprocess.CalledProcessError as e:
        print("\n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n  ")
        print(f"Standard output: \n {e.stdout} ")
        print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ")


############################################################################
#                                                                          #
#                         ▒▒▒▒▒▄██████████▄▒▒▒▒▒                           #
#                         ▒▒▒▄██████████████▄▒▒▒                           #
#                         ▒▒██████████████████▒▒                           #
#                         ▒▐███▀▀▀▀▀██▀▀▀▀▀███▌▒                           #
#                         ▒███▒▒▌■▐▒▒▒▒▌■▐▒▒███▒                           #
#                         ▒▐██▄▒▀▀▀▒▒▒▒▀▀▀▒▄██▌▒                           #
#                         ▒▒▀████▒▄▄▒▒▄▄▒████▀▒▒                           #
#                         ▒▒▐███▒▒▒▀▒▒▀▒▒▒███▌▒▒                           #
#                         ▒▒███▒▒▒▒▒▒▒▒▒▒▒▒███▒▒                           #
#                         ▒▒▒██▒▒▀▀▀▀▀▀▀▀▒▒██▒▒▒                           #
#                         ▒▒▒▐██▄▒▒▒▒▒▒▒▒▄██▌▒▒▒                           #
#                         ▒▒▒▒▀████████████▀▒▒▒▒                           #
#                          Program Begins Here                             #
#                                                                          #
############################################################################

command = './BlackHawk_tot.x parameter file'

error_list = []
current_dir = os.path.dirname(os.path.realpath(__file__))
spin_folder = 'Spin0.0'
output = '/spin00_mass'


main_folder_path = os.path.join(current_dir + '/src/tables/fM_tables/' + spin_folder )
destination_folder = os.path.join(current_dir +'/src/tables/fM_tables')
output_folder = current_dir + output
evolutions_folder = os.path.join(current_dir + '/results/test')
parameter_file = "parameters.txt"
fM_file_name = 'fM_add0.txt'
gM_file_name = 'gM_add0.txt'
# mdm_values = [0e+0, 1e-4,2e-4, 3e-4] 
mdm_values = [2e-4] 

for i in os.listdir(main_folder_path):
    if i != ".DS_Store":
        mass_folder =  os.path.join(output_folder + '/{:.1e}'.format(float(i)))
        os.makedirs(mass_folder, exist_ok=True)
        for j in os.listdir(main_folder_path + "/{}".format(i)):
            if j != ".DS_Store":
                dof_folder =  output_folder + '/{:.1e}'.format(float(i)) + '/{:.1e}'.format(float(j))
                os.makedirs(dof_folder, exist_ok=True)

                print("\n **************** {} ****************".format(('Mass','DoF')))
                print("**************** {} **************** \n ".format((i,j)))
                move_text_file(os.path.join(main_folder_path + "/{}".format(i)+ "/{}".format(j)), destination_folder, fM_file_name )
                move_text_file(os.path.join(main_folder_path + "/{}".format(i)+ "/{}".format(j)), destination_folder, gM_file_name )

                execute_BH(current_dir,parameter_file, i, j)

                move_text_file(destination_folder, os.path.join(main_folder_path + "/{}".format(i)+ "/{}".format(j)), fM_file_name )
                move_text_file(destination_folder, os.path.join(main_folder_path + "/{}".format(i)+ "/{}".format(j)), gM_file_name )
                print("\n ****************************************** \n\n")

                move_text_file(evolutions_folder, dof_folder, 'life_evolutions.txt' )

                old_filename =  dof_folder + "/life_evolutions.txt"
                new_filename =  dof_folder + "/Kerr_a0.8_{}M".format(i) +"_{}dof.txt".format(j)
                os.rename(old_filename, new_filename)
                remove_first_three_lines(new_filename)
print('Done')
for i in error_list:
    print('Error at {}'.format(i))
    


