import os

def generate_file_structure(directory, output_file):
    with open(output_file, 'w') as file:
        file.write('```mermaid\n')
        file.write('flowchart LR\n')  # Use flowchart instead of graph
        file.write('  linkStyle default interpolate basis\n')  # Default style
        file.write('  classDef default stroke-width:2px,stroke:#333\n')
        for root, dirs, files in os.walk(directory):
            dirs[:] = [d for d in dirs if not d.startswith('_')]
            parent = os.path.basename(root) or directory
            for d in dirs:
                if not d.startswith('_'):
                    file.write(f'    {parent} --> {d}\n')
            for f in files:
                file.write(f'    {parent} --> {f}\n')
        file.write('```\n')

if __name__ == "__main__":
    directory = "/home/esp/scriptDevelopment/drFrankenstein/Laboratory"
    output_file = "MODULE_STRUCTURE.md"
    generate_file_structure(directory, output_file)