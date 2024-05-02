error = function(e){ 
  message("Returned regression error")
  reg_error <<- e[1]
  failed <<- T
}