import os
import re
import json
import datetime
import pandas as pd
import SimpleITK as sitk
import argschema as ags


class IO_Schema(ags.ArgSchema):
    out_dir = ags.fields.OutputDir(required=True, description="Path to output directory", default= r'\\allen\programs\celltypes\workgroups\mousecelltypes\macaque_structures')


def process_json(jblob, annotation, structures ) :

    if jblob['markups'][0]['type'] != 'Fiducial':
        raise KeyError('Not a fiducial.')
    
    if 'controlPoints' not in jblob['markups'][0] :
        raise KeyError("No control points found.")
    
    locs = []
    for j in jblob['markups'][0]['controlPoints'] :

        info = {}
        info['specimen_name'] = re.sub(r'-(\d+)$', lambda m: f".{int(m.group(1)):02d}", j['label'])

        # Extract pin location 
        pos = j['position']
        info['x'] = pos[0]
        info['y'] = pos[1]
        info['z'] = pos[2]
        point = (info['x'], info['y'], info['z'])
        
        # -- this simply divides cooordinates by resolution/spacing to get the pixel index
        pixel = annotation.TransformPhysicalPointToIndex(point)
        sid = annotation.GetPixel(pixel)
        info['structure_id'] = sid
        
        if sid not in structures.index :
            print(info)
            print(f"WARNING: not a valid structure {sid}, skipping {info['specimen_name']}")
            continue
        
        info['structure_acronym'] = structures.loc[sid]['acronym']

        locs.append(info)

    return locs


def get_soma_and_fiducial_pins():

    #macaque annotation volume
    model_directory = r'\\allen\programs\celltypes\workgroups\humancelltypes\HMBA_annotations\macaque\annotations\Mac25Rhesus_v2'
    annotation_file = os.path.join(model_directory, 'Mac25Rhesus_v2.D99Atlas_v2_Subcortical_RegByRIKEN1_0.16mm.nii.gz')
    annotation = sitk.ReadImage( annotation_file )

    #macaque ccf structure data
    structure_file = os.path.join(model_directory, 'D99_Saleem_ITK_label_description.txt')
    structures = pd.read_csv(structure_file, sep=" ", header=None, names=['id', 'red', 'green', 'blue', 0, 1, 2, 'structure'])
    structures.set_index('id', inplace=True)
    structures[['acronym', 'name']] = structures['structure'].str.split(' - ', expand=True)

    #process all jsons (from all brains in the folder)
    cell_info = []
    json_root = r'\\allen\programs\celltypes\workgroups\mousecelltypes\Ingrid_Redford\NHP Pinning'
    for json_brain_dir, _, filenames in os.walk(json_root):
        for filename in filenames:
            if filename.endswith(".json"):
                json_path = os.path.join(json_brain_dir, filename)

                with open(json_path, 'r') as j_file:
                    jblob = json.load(j_file)

                processed = process_json(jblob, annotation, structures)
                cell_info.extend(processed)

    df = pd.DataFrame(cell_info)
    return df


def main(out_dir):
    print('\nAccumulating soma pins and finding associated structures...')

    df = get_soma_and_fiducial_pins(out_dir)

    print('Saving results...')
    x = datetime.datetime.now()
    fdt = str(x).split(' ')[0]
    y, m, d = fdt.split('-')
    dt = '{}{}{}'.format(y, m, d)

    output_file = "macaque_soma_structures.csv"
    df.to_csv(os.path.join(out_dir, output_file), index=False)

    output_file_archive = "macaque_soma_structures_{}.csv".format(dt)
    os.makedirs(os.path.join(out_dir, 'archive'), exist_ok=True)
    df.to_csv(os.path.join(out_dir, 'archive', output_file_archive), index=False)

    print('Done!')
    print(f'Results: {os.path.join(out_dir, output_file)}')

if __name__ == "__main__":
    parser = ags.ArgSchemaParser(schema_type=IO_Schema)
    main(parser.args['out_dir'])



