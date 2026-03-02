-- Adjust column widths for tbl-notation
function Table(tbl)
  if tbl.identifier == "tbl-notation" then
    -- Adjust column widths: narrower first column, wider second
    if tbl.colspecs and #tbl.colspecs == 3 then
      tbl.colspecs[1] = {pandoc.AlignDefault, 0.15}  -- Symbol column
      tbl.colspecs[2] = {pandoc.AlignDefault, 0.52}  -- Description column  
      tbl.colspecs[3] = {pandoc.AlignDefault, 0.28}  -- Defined in column
    end
  end
  return tbl
end
