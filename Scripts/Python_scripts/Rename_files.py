import os


def rename_files(directory):
    # Iterate over all files in the directory
    for filename in os.listdir(directory):
        # Construct full file path
        old_file = os.path.join(directory, filename)

        # Check if it is a file
        if os.path.isfile(old_file):
            # Replace '-' with '_' in the file name
            new_filename = filename.replace('-', '_')
            new_file = os.path.join(directory, new_filename)

            # Rename the file
            os.rename(old_file, new_file)
            print(f'Renamed: {old_file} -> {new_file}')


def rename_directories(directory):
    # Iterate over all items in the directory
    for item in os.listdir(directory):
        # Construct full item path
        old_item = os.path.join(directory, item)

        # Check if it is a directory
        if os.path.isdir(old_item):
            # Replace '-' with '_' in the directory name
            new_item_name = item.replace('-', '_')
            new_item = os.path.join(directory, new_item_name)

            # Rename the directory
            os.rename(old_item, new_item)
            print(f'Renamed: {old_item} -> {new_item}')


# Specify the directory you want to rename files in
directory_path = '/Users/at991/OneDrive - University of Exeter/Data/Cambridge_Project/CheckMbins/'
rename_files(directory_path)
