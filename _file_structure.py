import os

def generate_file_structure(directory, output_file):
    with open(output_file, 'w') as file:
        for root, dirs, files in os.walk(directory):
            level = root.replace(directory, '').count(os.sep)
            indent = ' ' * 4 * level
            file.write(f'{indent}- {os.path.basename(root)}/\n')
            sub_indent = ' ' * 4 * (level + 1)
            for f in files:
                file.write(f'{sub_indent}- {f}\n')

if __name__ == "__main__":
    directory = "/home/esp/scriptDevelopment/drFrankenstein/src"
    output_file = "module_structure.md"
    generate_file_structure(directory, output_file)