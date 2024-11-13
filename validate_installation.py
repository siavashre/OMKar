import filecmp
import os

def compare_directories_recursively(dir1, dir2):
    comparison = filecmp.dircmp(dir1, dir2)
    differences_found = False

    # Check for files and directories only in dir1 or dir2
    if comparison.left_only:
        differences_found = True
        print(f"File missing: ")
        for item in comparison.left_only:
            print(f"  {os.path.join(dir1, item)}")

    if comparison.right_only:
        differences_found = True
        print(f"File added: ")
        for item in comparison.right_only:
            print(f"  {os.path.join(dir2, item)}")

    ## do not check the log content
    if dir1.split('/')[-1] != 'logs':
        # Check for files that are present in both directories but differ in content
        for diff_file in comparison.diff_files:
            if diff_file.split('.')[-1] == 'pdf':
                continue
            differences_found = True
            print(f"Different content: {os.path.join(dir1, diff_file)} and {os.path.join(dir2, diff_file)}")

        # Check for files that are present in both directories and are identical
        for common_file in comparison.common_files:
            if common_file.split('.')[-1] == 'pdf':
                continue
            file1 = os.path.join(dir1, common_file)
            file2 = os.path.join(dir2, common_file)
            if not filecmp.cmp(file1, file2, shallow=False):  # Compare content fully
                differences_found = True
                print(f"Different content: {file1} and {file2}")

    # Recursively check subdirectories
    for sub_dir in comparison.common_dirs:
        sub_dir1 = os.path.join(dir1, sub_dir)
        sub_dir2 = os.path.join(dir2, sub_dir)
        if not compare_directories_recursively(sub_dir1, sub_dir2):
            differences_found = True

    return not differences_found

status = compare_directories_recursively('test_intended_output/', 'test_output/')
if status:
    print('Installation Validation: Passed')
else:
    print('Installation Validation: Failed')
