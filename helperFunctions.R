### Function to create directory names
makePath <- function(appendName=params$xtraName,name=""){
  xtraName <- appendName
  xtraName <- ifelse(xtraName=="",paste0(xtraName,"_"),paste0("_",xtraName,"_"))
  pathDir <- paste("results",ifelse(appendName=="","",paste0("_",appendName)),"/results",xtraName,name,collapse = "/",sep="")
  return(pathDir)
}

### String modifier functions
camel_case <- function(x,war=F) {
  if (war) cat ("*Forcing title to camelCase*\n")
  x <- gsub("(^|[^[:alnum:]])([[:alnum:]])", "\\U\\2", x, perl = TRUE)
  x <- paste(c(tolower(strsplit(x,"")[[1]][1]),strsplit(x,"")[[1]][-1]),collapse = "")
  return(x)
}
## Remove punctuations from string
clean_punct <- function(text){
  text <- gsub("[[:punct:]]", " ",text)
  return(text)
}
## Create a path to save a file
final_destination <- function(title,pathTo,prefix="",suffix="",extension="",doCamel=F){
  # Function to create the destination of a file, given a path (path/to/file), a title (titleOfFile) and an extension (.ext)
  # to form a path/to/file/titleOfFile.ext string to save the file.
  ## The function can accept a: 
  # prefix parameter (path/to/file/01prefix_titleOfFile.ext)
  # suffix parameter (path/to/file/titleOfFile_suffix.ext)
  # both (path/to/file/01prefix_titleOfFile_suffix.ext)
  ## 
  if (extension=="") warning("No extension provided\n")
  if (title=="") warning("No title provided\n")
  if (prefix !="") prefix <- paste0(prefix,"_")
  if (suffix !="") suffix <- paste0("_",suffix)
  if (extension !="") extension <- paste0(".",extension)
  ## Force to camel case in case there are spaces in the title
  if (doCamel) {title <- camel_case(clean_punct(title))}
  destPath <- paste(pathTo,"/",prefix,title,suffix,extension,sep="")
  return(destPath)
}
