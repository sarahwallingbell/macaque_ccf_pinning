@echo off
call C:\Users\sarah.wallingbell\AppData\Local\anaconda3\Scripts\activate.bat morph_utils_v4
python \\allen\programs\celltypes\workgroups\mousecelltypes\SarahWB\github_projects\macaque_ccf_pinning\macaque_ccf_pinning\macaque_soma_pin_structures.py
call conda deactivate