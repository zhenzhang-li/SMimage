.onLoad <- function(...){
	packageStartupMessage("\n")
	packageStartupMessage("Welcome to SMimage.")
	packageStartupMessage("\n")
	packageStartupMessage("Version: ",utils::packageDescription('SMimage')$Version)
	packageStartupMessage("\n")
	packageStartupMessage("If this is your first time running SMimage you should see ?SMimage")
	packageStartupMessage("\n")
	}