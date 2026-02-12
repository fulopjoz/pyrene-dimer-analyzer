! {header_line}

%pal
  nprocs {nprocs}
end

%maxcore {maxcore}

%tddft
  nroots {nroots}
end

%esd
  esdflag abs
  gshessian "{gshessian_file}"
  doht true
  lines gauss
  inlinew {inlinew}
  states {states}
end

* xyzfile {charge} {mult} {xyz_file}

