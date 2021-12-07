# DASiRe

# Serverside
Start the serverside app by running the following command from the top directory of this git.

`docker run -p 3838:3838 -v $(pwd)/serverside/R:/srv/shiny-server/ -d --rm dasire-server:0.1`
