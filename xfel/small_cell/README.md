notes for readme file

Find spots to prepare for powder pattern
With command.py containing:

```
import sys, os
i = int(sys.argv[1])
j = int(sys.argv[2])
cmd="mpirun -n 64 cctbx.small_cell_process spotfind.phil /net/dials/raid1/elyseschriber/sacla2019/data/data/78%d-%d/run78%d-%d.h5 output.output_dir=%s mp.method=mpi"
odir="%d-%d"%(i,j)
cmd_ji= cmd %(i,j, i, j, odir)
print (cmd_ji)
os.mkdir (odir)
os.system (cmd_ji)
```

and spotfind.phil containing:
```
input.reference_geometry=/net/dials/raid1/elyseschriber/processing/protk/dan/step2/step2_refined2.expt
dispatch.hit_finder {
  minimum_number_of_reflections = 3
  maximum_number_of_reflections = 40
}
dispatch.index=False
output.composite_output=True
spotfinder.filter.min_spot_size=3
```

do:
```
$ for n in {3192..3200}; do for nn in {0..2}; do python command.py $n $nn; done; done
```

Then combine experiments:
```
$ mkdir combined
$ for d in 3*; do dials.combine_experiments $d/*imported.expt $d/*strong.refl output.experiments=combined/$d.expt output.reflections=combined/$d.refl reference_from_experiment.detector=0 & done
$ cd combined; dials.combine_experiments *.expt *.refl reference_from_experiment.detector=0
```

Make powder pattern:
```
cctbx.xfel.powder_from_spots combined.* output.xy_file=powder.xy
```

Manually pick peaks and enter them in peaks.txt, then do:
```
cctbx.xfel.candidate_cells nproc=64 input.peak_list=peaks_new.txt input.powder_pattern=powder.xy search.timeout=300
```

Take this indexing phil file and make 5 copies:

```
input.reference_geometry=/net/dials/raid1/elyseschriber/processing/protk/dan/step2/step2_refined2.expt
dispatch {
  hit_finder {
    minimum_number_of_reflections = 3
    maximum_number_of_reflections = 40
   }
  }
dispatch {
  refine = True
}
output {
  composite_output=True
  experiments_filename = None
  strong_filename = None
}

spotfinder {
  filter {
    min_spot_size = 3
}
refinement {
  parameterisation {
    crystal {
      fix = all *cell orientation
    }
  }
  reflections {
    outlier {
      algorithm = null auto mcd tukey sauter_poon
    }
  }
}
small_cell {
  powdercell = "1,2,3,90,90,90"
  spacegroup = "P2"
  high_res_limit = 1.2
  min_spots_to_integrate = 3
  faked_mosaicity = 0.1
  spot_connection_epsilon = 0.01
  d_ring_overlap_limit = None
}
```

Name the copies 1.phil through 5.phil. Manually enter the unit cells and space groups of the top 5 candidates from candidate_cells.

With run.py containing:
```
import sys, os
phil_root = sys.argv[1]
i, j = 3192, 0
cmd="mpirun -n 64 cctbx.small_cell_process %s.phil /net/dials/raid1/elyseschriber/sacla2019/data/data/78%d-%d/run78%d-%d.h5 output.output_dir=%s mp.method=mpi"
cmd_ji= cmd %(phil_root, i,j, i, j, phil_root)
print (cmd_ji)
os.mkdir (phil_root)
os.system (cmd_ji)
```

do

```
$ for n in 1 2 3 4 5; do python run.py $n; done
```
and then
```
$ for n in {1..5}; do echo;  echo $n; egrep 'powdercell|spacegroup' $n.phil;  grep real_space_a $n/*refined.expt 2>/dev/null|wc -l; done
```


After choosing the correct cell, prepare this file integrate.phil:
```
input.reference_geometry=/net/dials/raid1/elyseschriber/processing/protk/dan/step2/step2_refined2.expt
dispatch {
  hit_finder {
    minimum_number_of_reflections = 3
    maximum_number_of_reflections = 40
   }
  }
dispatch {
  refine = True
}
output {
  composite_output=True
  experiments_filename = None
  strong_filename = None
}

spotfinder {
  filter {
    min_spot_size = 3
  }
}
refinement {
  parameterisation {
    crystal {
      fix = all *cell orientation
    }
  }
  reflections {
    outlier {
      algorithm = null auto mcd tukey sauter_poon
    }
  }
}
small_cell {
  powdercell = "5.93762, 7.32461, 29.202, 90, 95.4404, 90"
  spacegroup = "C2"
  high_res_limit = 0.8
  min_spots_to_integrate = 3
  faked_mosaicity = 0.1
  spot_connection_epsilon = 0.005
  d_ring_overlap_limit = None
}
```

And this script run.py:
```
import sys, os
phil_root = sys.argv[1]
i = int(sys.argv[2])
j = int(sys.argv[3])
output_dir = "{}-{}-{}".format(phil_root, i, j)
logging_dir = output_dir + '/stdout'
cmd = "mpirun -n 60 cctbx.small_cell_process {}.phil /net/dials/raid1/elyseschriber/sacla2019/data/data/78{}-{}/run78{}-{}.h5 output.output_dir={} output.logging_dir={} mp.method=mpi"
cmd_ji = cmd.format(phil_root, i, j, i, j, output_dir, logging_dir)
print (cmd_ji)
os.mkdir (output_dir)
os.mkdir(logging_dir)
os.system (cmd_ji)
```
Prepare a todo script and run it: 
```
$ for m in {3133..3214}; do for n in {0..2}; do echo "python run.py integrate $m $n" >> todo.sh; done; done
$ source todo.sh
```
After integration, merging is performed in two steps:

```
$ mkdir merging; cd merging
$ cat > mark1.phil
input.path=../3*
dispatch.step_list = input balance model_scaling modify errors_premerge scale statistics_unitcell statistics_beam model_statistics statistics_resolution group errors_merge statistics_intensity merge statistics_intensity_cxi
scaling.algorithm=mark1
scaling.unit_cell=5.93762 7.32461 29.202 90 95.4404 90
scaling.space_group=C2
merging.d_min=1.15
merging.merge_anomalous=True
merging.error.model=errors_from_sample_residuals
statistics.n_bins=20
output.prefix=mark1
output.do_timing=True
output.output_dir=.

$ cat > mark0.phil
input.path=../3*
dispatch.step_list = input balance model_scaling modify errors_premerge scale statistics_unitcell statistics_beam model_statistics statistics_resolution group errors_merge statistics_intensity merge statistics_intensity_cxi
scaling.model=mark1_all.mtz
merging.d_min=1.15
merging.merge_anomalous=True
output.do_timing=True
statistics.n_bins=20
merging.error.model=errors_from_sample_residuals
output.prefix=mark0
output.do_timing=True
output.output_dir=.

$ mpirun -n 32 cctbx.xfel.merge mark1.phil
$ mpirun -n 32 cctbx.xfel.merge mark0.phil
$ iotbx.reflection_file_converter mark0_all.mtz --label="Iobs,SIGIobs" --shelx=mark0.hkl
```

At this point, the file mark0.hkl can be used for structure solution. To perform a final round of scaling with a partially completed reference structure (here combined with reindexing to resolve an indexing ambiguity), save the structure in standard cif format and do this:
```
$ cat > rtr.phil
input.path=../3*
dispatch.step_list = input balance model_scaling modify modify_reindex_to_reference errors_premerge scale statistics_unitcell statistics_beam model_statistics statistics_resolution group errors_merge statistics_intensity merge statistics_intensity_cxi
scaling.model=full.cif
merging.d_min=1.15
merging.merge_anomalous=True
output.do_timing=True
statistics.n_bins=20
merging.error.model=errors_from_sample_residuals
output.prefix=rtr
output.do_timing=True
output.output_dir=.

$ mpirun -n 32 cctbx.xfel.merge rtr.phil
$ iotbx.reflection_file_converter rtr_all.mtz --label="Iobs,SIGIobs" --shelx=rtr.hkl
```

