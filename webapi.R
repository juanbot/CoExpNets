
globalReportOnGenesWEB = function(){
  req <- curl::curl_fetch_memory(paste0("https://snca.atica.um.es/api_test/globalReportOnGenes?tissues=",
    paste0(tissues,collapse=","),
    "&categories=",paste0(categories,collapse=","),
    "&genes=",paste0(genes,collapse=","),"/"))

  myjson = jsonlite::fromJSON(rawToChar(req$content))
  return(myjson)
}
