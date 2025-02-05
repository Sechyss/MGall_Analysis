import xml.etree.ElementTree as ET

def read_xml(file_path):
    tree = ET.parse(file_path)
    root = tree.getroot()
    
    variables = {}
    for child in root:
        variables[child.tag] = child.text
    
    return variables

if __name__ == "__main__":
    file_path = 'path_to_your_xml_file.xml'
    variables = read_xml(file_path)
    print(variables)