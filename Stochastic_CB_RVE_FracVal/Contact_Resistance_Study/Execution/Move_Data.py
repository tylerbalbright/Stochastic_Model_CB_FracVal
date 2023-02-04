#!/usr/bin/env python3
import os
import csv
import shutil

temp_dirs = []
perm_dirs = []

cwd = os.getcwd()
source_dir = (cwd + "/RVE_Data/")
for temp_folders in os.listdir(source_dir):
    if temp_folders[0] == "_": continue
    basename, decimal = os.path.splitext(temp_folders)
    if temp_folders[0] == ".": continue
    #XY, Wcb = basename.split('_')
    temp_dirs.append(basename+decimal)

storage_dir = "/state/partition4/Tyler/Conductivity_Study/Strain_Study/RVE_Data/"
for perm_folders in os.listdir(storage_dir):
    if perm_folders[0] == "_": continue
    basename, decimal = os.path.splitext(perm_folders)
    if perm_folders[0] == ".": continue
    #XY, Wcb = basename.split('_')
    perm_dirs.append(basename+decimal)

perm_file_types = []
perm_file_nums = []
dir_exists = False
for i in range(len(temp_dirs)):
    dir_exists = False
    for j in range(len(perm_dirs)):
        if temp_dirs[i] == perm_dirs[j]:
            # the folders match, add temp files to the storage bin
            # by checking to see if similar files exist!
            dir_exists = True
            for filename in os.listdir(storage_dir + perm_dirs[j]):
                if filename[0] == "_": continue
                basename1, extension1 = os.path.splitext(filename)
                basename, extension = os.path.splitext(basename1)
                thickness, DQ, Fill, num = basename.split('_')
                temp_type = thickness + "_" + DQ + "_" + Fill + "_"
                if temp_type in perm_file_types:
                    perm_file_nums[perm_file_types.index(temp_type)] += 1
                else:
                    perm_file_types.append(temp_type)
                    perm_file_nums.append(0)

            for filename in os.listdir(source_dir + temp_dirs[i]):
                # Loop through the temp directory and move the files
                # to the permanent directory.
                if filename[0] == "_": continue
                basename1, extension1 = os.path.splitext(filename)
                basename, extension = os.path.splitext(basename1)
                try:
                    ignore, thickness, DQ, Fill, num = basename.split('_')
                except:
                    print(basename)
                    exit()
                temp_type = thickness + "_" + DQ + "_" + Fill + "_"
                if temp_type in perm_file_types:
                    perm_file_nums[perm_file_types.index(temp_type)] += 1
                    old_name = source_dir + temp_dirs[i] + "/" + filename
                    new_name = storage_dir + perm_dirs[j] + "/" + temp_type + str(perm_file_nums[perm_file_types.index(temp_type)]) + ".xyz.xz"
                    shutil.move(old_name, new_name)
                else:
                    old_name = source_dir + temp_dirs[i] + "/" + filename
                    filename_edit = filename[len(ignore)+1:]
                    new_name = storage_dir + perm_dirs[j] + "/" + filename_edit
                    perm_file_types.append(temp_type)
                    perm_file_nums.append(0)
                    shutil.move(old_name ,new_name)
    if dir_exists is False:
        os.mkdir(storage_dir + temp_dirs[i])
        moved_files = []
        for filename in os.listdir(source_dir + temp_dirs[i]):
            # Loop through the temp directory and move the files
            # to the permanent directory.
            basename1, extension1 = os.path.splitext(filename)
            basename, extension = os.path.splitext(basename1)
            ignore, thickness, DQ, Fill, num = basename.split('_')
            num_int = 0
            temp_name = thickness + "_" + DQ + "_" + Fill + "_" + str(num_int) + ".xyz.xz"
            for z in range(len(moved_files)):
                if temp_name == moved_files[z]:
                    num_int = num_int + 1
                    temp_name = thickness + "_" + DQ + "_" + Fill + "_" + str(num_int) + ".xyz.xz"
            temp_name = thickness + "_" + DQ + "_" + Fill + "_" + str(num_int) + ".xyz.xz"
            print(temp_name)
            moved_files.append(temp_name)
            old_name = source_dir + temp_dirs[i] + "/" + filename
            new_name = storage_dir + temp_dirs[i] + "/" + temp_name
            shutil.move(old_name, new_name)
    os.rmdir(source_dir + temp_dirs[i])
