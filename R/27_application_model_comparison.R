ROOT <- file.path(RESULT_DIR,'weighted_single_sample')
for(p in c('tables','workflow','figures'))dir.create(file.path(ROOT,p),recursive=TRUE,showWarnings=FALSE)
source(file.path(SOURCE_ROOT,'application_models','prepare_data.R'),local=environment())
source(file.path(SOURCE_ROOT,'application_models','compare_models.R'),local=environment())
