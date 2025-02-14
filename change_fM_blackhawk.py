import subprocess
import shutil
import time
import os
import re
import fileinput
import numpy as np 

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

# def update_spin_dm(fmc_path, new_value):
#     with open(fmc_path, "r") as file:
#         content = file.read()

#     # Regex pattern to match 'double spin_DM = any_value;'
#     pattern = r"(double spin_DM\s*=\s*)[-+]?[0-9]*\.?[0-9]+;"

#     # Correct replacement using `.format()`
#     updated_content = re.sub(pattern, r"\1{};".format(new_value), content)

#     # Write the updated content back to the file
#     with open(fmc_path, "w") as file:
#         file.write(updated_content)


def execute_fM(fM_name,gM_name,mass,dof):
    try:
        # Step 1: Run 'make' to compile fM.c
        subprocess.run(['make'], cwd=fM_folder, check=True)
        print("\n Compilation successful. \n")

        # Step 2: Run './fM.x 0' and handle the input prompts
        process = subprocess.Popen(['./fM.x', '0'], cwd=fM_folder, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


        # Send input for the two prompts
        input_data = f"{fM_name}\n{gM_name}\n"
        output, error = process.communicate(input=input_data)

        # Check for errors during execution
        if process.returncode != 0:
            print(f"Error occurred: {error}")
            error_list.append((mass,dof))
        else:
            print("Execution completed successfully.")
            print("Output:")
            print(output)

    except subprocess.CalledProcessError as e:
        print(f"An error occurred during the process: {e}")


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

################### INPUT VARIABLES ###################

mdm_values = [0e+0, 1e-4,2e-4, 3e-4] 
# mdm_values = [2e-4] 
# mdm_values = [3e-4,4e-4,5e-5] 
# dof_values = [1e01,1e2,1e3,1e4,1e05,1e7,1e9,1e11,1e12,1e15,1e20,1e25,1e30,1e40,1e50]
# dof_values = [1e01,1e2,1e3,1e4,1e05,1e7,1e9]
dof_values = [1e0,1e1,1e2,1e3,1e4,1e5]

dm_spin = 0
spin_folder = 'Spin0.0'


#######################################################



######################## PATHS ########################

current_dir = os.path.dirname(os.path.realpath(__file__))
fM_file_name = 'fM_add0.txt'
gM_file_name = 'gM_add0.txt'
fM_folder = os.path.join(current_dir,"scripts/greybody_scripts/fM")

fM_path = fM_folder + fM_file_name
gM_path = fM_folder + gM_file_name
fmc_path = fM_folder + "/fM.c"

#######################################################

# update_spin_dm(fmc_path, dm_spin)

# Change spin to what we want in fM.c
error_list = []


## FIX FOLDER PATH

for i in mdm_values:

    destination_folder = os.path.join(current_dir,'src/tables/fM_tables/'+ spin_folder +'/{:.0e}'.format(float(i)))
    os.makedirs(destination_folder, exist_ok=True)

    with open(fmc_path, "r") as file:
        lines = file.readlines()

# Modify the line that contains 'double m_DM'
    with open(fmc_path, "w") as file:
        for line in lines:
            if re.match(r"^\s*double\s+m_DM\s*=", line):  # Match 'double m_DM = ...'
                file.write(f"double m_DM = {i};\n")  # Replace with new value
            else:
                file.write(line)

    print("Updated `m_DM = {}` successfully!".format(i))
    
    for j in dof_values:
        mass_folder =  os.path.join(destination_folder,'{:.1e}'.format(float(j)))
        if os.path.exists(mass_folder) and os.path.isdir(mass_folder):
            print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n ")
            print(f"Folder '{mass_folder}' exists. Skipping process. \n")
            print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        else:
            print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n ")
            print(f"Folder '{mass_folder}' does not exist. Running process...")

           
            os.makedirs(mass_folder, exist_ok=True)
            with open(fmc_path, "r") as file:
                lines = file.readlines()

        # Modify the line that contains 'double dof_DM'
            with open(fmc_path, "w") as file:
                for line in lines:
                    if re.match(r"^\s*double\s+dof_DM\s*=", line):  # Match 'double m_DM = ...'
                        file.write(f"double dof_DM = {j};\n")  # Replace with new value
                    else:
                        file.write(line)

            print("Updated `dof_DM = {}` successfully! \n".format(j))

            #Create both the fM,gM folders 
            execute_fM(fM_file_name,gM_file_name,i,j)

            # Move both fM,gM
            move_text_file(fM_folder, mass_folder,fM_file_name)
            move_text_file(fM_folder, mass_folder,gM_file_name)        
            print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Finished with M_dm: {}, ".format(i) + " DoF: {} %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%".format((j)))

            

print('Number of errors = {} \n'.format(len(error_list)))
if len(error_list) > 0:
    for er in error_list:
        print(er)


# pipeline will be to change mdM first, then a series of DoF in fM.c +
# create fM_table +
# move fM_table to new folder
# run this for 