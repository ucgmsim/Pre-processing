stat_dir = os.path.join(global_root, 'StationInfo')
stat_file = os.path.join(stat_dir, 'cantstations.ll')
stat_coords = os.path.join(stat_dir, 'fd_%s.statcords'%sufx)
FD_STATLIST = os.path.join(stat_dir,'fd_%s.ll'%sufx)
stat_vs_ref = os.path.join(stat_dir,'cantstations_cant1D_v1.vs30ref')# station file with vs30 reference values
stat_vs_est = os.path.join(stat_dir,'cantstations.vs30')# station file with vs30 estimates

