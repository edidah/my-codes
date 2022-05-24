install.packages("rbibutils")
library(rbibutils)

##single file conversion 

bibConvert("Bloom et al., 2021.nbib",  "Bloom et al., 2021.ris")


#multiple file conversion steps
## 1. subset the files with the specific extension and list them
test_files <- list.files(
  path = "set file path", # replace with the directory you want
  pattern = ".*nbib", # has ".nbib",
  full.names = TRUE # include the directory in the result
)

##Optional; Determining the number of files
length(test_files)

### optional; if renaming otherwise proceed to next step
no_file_ext <- tools::file_path_sans_ext(test_files)
no_file_ext
new_names <- sprintf("%s.nbib", no_file_ext, seq_along(no_file_ext))
new_names[1]

###convert using bibConvert function and iterate through the list to generate the converted list
for (i in seq_along(test_files)) {
  bibConvert(test_files[[i]], paste0("test_file", i, ".ris")) #iterate and write file
}


