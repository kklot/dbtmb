#' Roxygen commands
#'
#' @useDynLib dbtmb
#'
dummy <- function(){
  return(NULL)
}

# https://stackoverflow.com/a/43069470
code <- "int wstat; while (waitpid(-1, &wstat, WNOHANG) > 0) {};"
#'
#' @export
wait_children_kill <- inline::cfunction(
  body = code, includes = "#include <sys/wait.h>", convention = ".C"
)
