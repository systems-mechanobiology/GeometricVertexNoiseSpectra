-- Lua filter to handle acronyms for HTML output
-- Converts \ac{ACRONYM} to the acronym text for HTML
-- Leaves LaTeX output unchanged (handled by acro package)

local acronyms = {
  SC = "SC",
  SSI = "SSI",
  S3I = "S3I",
  TT = "TT",
  CTPL = "CTPL"
}

local acronyms_long = {
  SC = "stochastic co-spectrality",
  SSI = "stochastic spectral identifiability",
  S3I = "stochastic spectral separation index",
  TT = "tempered tail",
  CTPL = "calibrated tempered power-law"
}

local first_use = {}

function RawInline(el)
  if el.format == "tex" and el.text:match("^\\ac{") then
    -- Extract acronym name
    local acr = el.text:match("\\ac{([^}]+)}")
    if acr and acronyms[acr] then
      -- For non-LaTeX formats (HTML), expand the acronym
      if FORMAT:match("html") then
        if not first_use[acr] then
          first_use[acr] = true
          return pandoc.Str(acronyms_long[acr] .. " (" .. acronyms[acr] .. ")")
        else
          return pandoc.Str(acronyms[acr])
        end
      end
    end
  end
  return el
end
