-- strip-linebreaks.lua
-- Converts hard line breaks (from trailing double spaces in Markdown)
-- to soft breaks, preventing unwanted line breaks in PDF output.
-- One-sentence-per-line source style is preserved; Pandoc joins
-- consecutive non-blank lines with a space as usual.

function LineBreak()
  return pandoc.SoftBreak()
end
