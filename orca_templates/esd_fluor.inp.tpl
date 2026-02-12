! {header_line}

%pal
  nprocs {nprocs}
end

%maxcore {maxcore}

%tddft
  nroots {nroots}
  iroot {iroot}
end

%esd
  esdflag fluor
  gshessian "{gshessian_file}"
  eshessian "{eshessian_file}"
  doht true
  lines voigt
  linew {linew}
  inlinew {inlinew}
  states {states}
end

* xyzfile {charge} {mult} {xyz_file}

