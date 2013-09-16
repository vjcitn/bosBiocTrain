laconv = function(x) {
 x = gsub("\\\\section\\{(.*)\\}", "<h2>\\1</h2>", x)
 x = gsub("\\\\texttt\\{(.*)\\}", "<code>\\1</code>", x)
 x = gsub("\\\\verb+(.*)+\\}", "<code>\\1</code>", x)
 x = gsub("\\\\emph\\{(.*)\\}", "<emph>\\1</emph>", x)
 x = gsub("\\\\subsection\\{(.*)\\}", "<h3>\\1</h3>", x)
 x = gsub("\\\\subsubsection\\{(.*)\\}", "<h4>\\1</h4>", x)
 x = gsub("\\\\begin\\{itemize\\}", "<ul>", x)
 x = gsub("\\\\begin\\{verbatim\\}", "<code>", x)
 x = gsub("\\\\end\\{verbatim\\}", "</code>", x)
 x = gsub("\\\\end\\{itemize\\}", "</ul>", x)
 x = gsub("\\\\item ", "<li>", x)
 x = gsub("<<(.*)>>=", "```{r \\1}", x)
 x = gsub("^@$", "```", x)
 x
}

convla = function(fn) {
  bk = readLines(fn)
  writeLines(bk, con=paste0(fn, ".bak"))
  writeLines(laconv(bk), con=paste0(fn, ".Rmd"))
}
