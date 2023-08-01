#-----------------------------Generic functions--------------------------------#

### Function to ensure directory exists before creating files in it 
ensure_directory <- function(directory){
  if(!dir.exists(directory)){
    dir.create(directory);
  }
}

### Function to create logger function 
make_logger <- function(filename, sep="\n"){
  if(file.exists(filename)){
    file.remove(filename);
  }
  function(...){
    text <- sprintf(...);
    cat(text, file=filename, sep=sep, append=T);
  }
}
