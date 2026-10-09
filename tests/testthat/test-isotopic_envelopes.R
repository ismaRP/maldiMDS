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

test_that("no mods returns formula", {
  form = peptide_formula('GPPGPPGPP')
  modform = add_modifications(form, NULL)

  expect_equal(form, modform)

})

test_that("correct simple formula", {
  form = peptide_formula('GQP')
  expected_form = 'C12H20N4O5'
  expect_equal(form, expected_form)

  mods = list(hyd=1, deam=1)
  modform = add_modifications(form, mods)
  expected_modform = 'C12H20N3O7'
  expect_equal(modform, expected_modform)

})


test_that("correct peptide formula", {
  form = peptide_formula('GVQGPPGPAGPR')
  expected_form = 'C47H76N16O14'
  expect_equal(form, expected_form)

  mods = list(hyd=1, deam=1)
  modform = add_modifications(form, mods)

  expected_modform = 'C47H76N15O16'
  expect_equal(modform, expected_modform)
})


test_that("correct isotopic isotopic substitutions", {
  iso_variants = get_isotopic_variants('C3H1O1S1', threshold=1e-3)

  agg_iso = get_n_isosubs(iso_variants)

  expect_length(agg_iso, 17)
  expect_equal(
    agg_iso,
    c(0,1,1,1,1,2,2,2,2,2,3,3,3,3,4,4,4)
  )
})






