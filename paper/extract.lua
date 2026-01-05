function Pandoc(doc)
  -- This creates a new bibliography containing only used references
  local refs = pandoc.utils.references(doc)
  -- We return a "fake" document that just contains those references
  return pandoc.Pandoc({}, {references = refs})
end
