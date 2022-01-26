####start the project
library("workflowr")

wflow_status()
############
wflow_git_config(user.name = "baigal628", user.email = "gali.bai@hotmail.com")
wflow_use_github("baigal628")
#################
wflow_build()
wflow_publish(c("analysis/multi-scatac_stepbystep.Rmd*", "docs/assets/*"),
              "Fix samples.json help information format")

wflow_publish("docs/assets/*",
              "Fix samples.json help information format")

wflow_git_push(dry_run = TRUE)
wflow_git_push()

### Add new file
wflow_open("analysis/snrna.Rmd")
