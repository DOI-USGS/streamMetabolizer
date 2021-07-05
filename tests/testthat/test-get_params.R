context('get_params')

dat <- data_metab('3','15')
mm <- metab_mle(specs("m_np_oi_tr_plrckm.nlm"), data=dat)
dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)), K600.daily=21)
mm2 <- metab_mle(specs("m_np_oi_tr_plrckm.nlm"), data=dat, data_daily=dat_daily)

test_that('get_params options are honored (for MLE models): uncertainty', {
  # uncertainty
  ps <- get_params(mm, uncertainty='sd')
  expect_equal(grep('\\.sd$', names(ps), value=TRUE), c('GPP.daily.sd','ER.daily.sd','K600.daily.sd'))
  ps <- get_params(mm, uncertainty='ci')
  expect_equal(grep('\\.lower$', names(ps), value=TRUE), c('GPP.daily.lower','ER.daily.lower','K600.daily.lower'))
  ps <- get_params(mm, uncertainty='none')
  expect_equal(length(grep('lower$|upper$|sd$', names(ps), value=TRUE)), 0)
})

test_that('get_params options are honored (for MLE models): fixed', {
  # fixed
  ps <- get_params(mm2, fixed='none')
  expect_equal(length(grep('\\.fixed$', names(ps), value=TRUE)), 0)
  expect_equal(length(grep('\\*', ps)), 0)
  ps <- get_params(mm2, fixed='columns')
  expect_equal(grep('\\.fixed$', names(ps), value=TRUE), c('GPP.daily.fixed','ER.daily.fixed','K600.daily.fixed'))
  ps <- get_params(mm2, fixed='stars')
  expect_true(all(sapply(dplyr::select(ps, -date), is.character)))
  expect_equal(grep('\\*', ps), match('K600.daily', names(ps)))
})

test_that('get_params options are honored (for MLE models): uncertainty+fixed', {
  # uncertainty + fixed
  ps <- get_params(mm, uncertainty='ci', fixed='stars')
  expect_equal(grep('\\.lower$', names(ps), value=TRUE), c('GPP.daily.lower','ER.daily.lower','K600.daily.lower'))
  expect_equal(grep('\\.fixed$', names(ps), value=TRUE), character(0))
  expect_equal(length(grep('\\*', ps)), 0)
  ps <- get_params(mm, uncertainty='ci', fixed='columns')
  expect_equal(grep('\\.lower$', names(ps), value=TRUE), c('GPP.daily.lower','ER.daily.lower','K600.daily.lower'))
  expect_equal(grep('\\.fixed$', names(ps), value=TRUE), c('GPP.daily.fixed','ER.daily.fixed','K600.daily.fixed'))
  ps <- get_params(mm2, uncertainty='ci', fixed='stars')
  expect_equal(grep('\\.lower$', names(ps), value=TRUE), c('GPP.daily.lower','ER.daily.lower'))
  expect_equal(grep('\\*', ps), match('K600.daily', names(ps)))
  ps <- get_params(mm2, uncertainty='ci', fixed='columns')
  expect_equal(grep('\\.lower$', names(ps), value=TRUE), c('GPP.daily.lower','ER.daily.lower'))
  expect_equal(grep('\\.fixed$', names(ps), value=TRUE), c('GPP.daily.fixed','ER.daily.fixed','K600.daily.fixed'))
})

test_that('get_params options are honored (for MLE models): messages', {
  # messages
  expect_true(all(c('warnings','errors') %in% names(get_params(mm, messages=TRUE))))
  expect_true(!any(c('warnings','errors') %in% names(get_params(mm, messages=FALSE))))
})

test_that('get_params works for each model type, basic GPP & ER equations', {
  dat <- data_metab('1','30')

  # empty model
  mm <- metab_model()
  expect_null(get_params(mm))

  # metab_mle
  mm <- metab_mle(data=dat)
  ps <- get_params(mm, uncertainty='none', messages=FALSE)
  expect_equal(names(ps), c('date','GPP.daily','ER.daily','K600.daily'))

  # metab_night
  mm <- metab_night(data=dat)
  ps <- get_params(mm, uncertainty='none', messages=FALSE)
  expect_equal(names(ps), c('date','ER.daily','K600.daily'))

  # metab_bayes
  mm <- metab_bayes(specs("b_np_oi_tr_plrckm.stan", burnin_steps=50, saved_steps=50, n_chains=1, n_cores=1), data=dat)
  pn <- get_param_names(mm)
  expect_equal(pn$required, c('GPP.daily','ER.daily','K600.daily'))
  ps <- get_params(mm, uncertainty='none', messages=FALSE)
  expect_equal(names(ps), c('date','GPP.daily','ER.daily','K600.daily'))

  # metab_sim
  dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)), GPP.daily=2, ER.daily=-3, K600.daily=21)
  mm <- metab_sim(specs(mm_name('sim'), K600_daily=NULL, GPP_daily=NULL, ER_daily=NULL), data=dat, data_daily=dat_daily)
  ps <- get_params(mm)
  expect_equal(names(ps), c('date','K600.daily','GPP.daily','ER.daily','err.obs.sigma','err.obs.phi','err.proc.sigma','err.proc.phi','discharge.daily'))
  ps <- get_params(mm, fixed='stars')
  expect_equal(grep('\\*', ps), match(c('K600.daily','GPP.daily','ER.daily'), names(ps)))
  dat_daily <- data.frame(date=as.Date(paste0("2012-09-", 18:20)), Pmax=6, alpha=0.001, ER.daily=-3, K600.daily=21)
  mm <- metab_sim(specs(mm_name('sim', GPP_fun='satlight'), Pmax=NULL, alpha=NULL, K600_daily=NULL, ER_daily=NULL), data=dat, data_daily=dat_daily)
  ps <- get_params(mm, fixed='columns')
  expect_equal(names(ps)[1:9], c('date','K600.daily','K600.daily.fixed','Pmax','Pmax.fixed','alpha','alpha.fixed','ER.daily','ER.daily.fixed'))
})
