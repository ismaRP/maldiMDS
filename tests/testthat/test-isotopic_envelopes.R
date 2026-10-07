test_that("long format isotopic envelopes have right shape", {
  seqs = c('GPPGPQGPP', 'GQPGPQGPP')
  ndeam = 2
  nhyd = c(2, 1)
  iso_env = get_isodists(seqs, ndeam, nhyd, sum, T)
  expect_equal(nrow(iso_env), 10)
  expect_equal(ncol(iso_env), 4 + ndeam)

  pep_idx = rep(c(1,2), each=5)
  mass_pos = rep(c(1,2,3,4,5), times=2)

  expect_equal(iso_env$pep_idx, pep_idx)
  expect_equal(iso_env$mass_pos, mass_pos)
})

test_that('correct envelopes', {
  test_envelope = readRDS(test_path('fixtures', 'test_envelope.rds'))
  iso_env = get_isodists(
    test_envelope$seqs,
    test_envelope$ndeam,
    test_envelope$nhyd,
    sum, T)
  expected_env = test_envelope$iso_env

  expect_equal(
    iso_env,
    expected_env
  )

})

test_that("deamidated isotopic envelopes are shifted", {
  seqs = c('GPPGPQGPP', 'GQPGPQGPP')
  ndeam = 2
  nhyd = c(2, 1)
  iso_env = get_isodists(seqs, ndeam, nhyd, sum, T)

  expect_equal(iso_env$deam_1[1], 0)

  expect_all_equal(iso_env$deam_2[1:5], 0)
  expect_all_equal(iso_env$deam_2[6:7], 0)

})
