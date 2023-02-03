Traceback (most recent call last):
  File "/scratch/bweight/software/SQD/src//NAMD/NAMD_main.py", line 129, in <module>
    main()
  File "/scratch/bweight/software/SQD/src//NAMD/NAMD_main.py", line 85, in main
    output.save_data(DYN_PROPERTIES)
  File "/gpfs/fs2/scratch/bweight/software/SQD/src/NAMD/output.py", line 73, in save_data
    if ( DYN_PROPERTIES['NAMD_METHOD'] in ["GFSH"] ):
KeyError: 'ACTIVE_STATE'
