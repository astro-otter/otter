# CHANGELOG

<!-- version list -->

## v1.1.1 (2026-02-19)

### Bug Fixes

- Upperlimit flagging issue
  ([`c91a709`](https://github.com/astro-otter/otter/commit/c91a709a780586ec157ae5470a164b92e90a1b5b))


## v1.1.0 (2026-01-27)

### Bug Fixes

- A NaN bug in _merge_arbitrary
  ([`b74436e`](https://github.com/astro-otter/otter/commit/b74436e00c1274524e89a69344a98025c312d20b))

- Edit pandas version requirement for wiserep_api
  ([`52ad10f`](https://github.com/astro-otter/otter/commit/52ad10f7364176548222bc56fe2695152de5d17d))

- Minor issue with pytest decorator
  ([`1b1dbe2`](https://github.com/astro-otter/otter/commit/1b1dbe23f56e66b04f452364ed6c2a39b840d2d8))

- Raising the query wiserep exception to keep track of the full traceback
  ([`4259496`](https://github.com/astro-otter/otter/commit/4259496e04b625d5f58c3b1746bc3e0988941b6e))

### Features

- Acatch an re-raise a better error if the object is not on wiserep
  ([`2036c21`](https://github.com/astro-otter/otter/commit/2036c2119f6c9cd0c5a3c0e942c09deb1518112c))


## v1.0.2 (2026-01-19)

### Bug Fixes

- Catch and handle documents that don't have the telescope column when cleaning the phot
  ([`11a9475`](https://github.com/astro-otter/otter/commit/11a947532df5e251f4b96ae4727790a9ccdb4bde))


## v1.0.1 (2026-01-19)

### Bug Fixes

- A bug in the photometry deduplication code
  ([`0227fb0`](https://github.com/astro-otter/otter/commit/0227fb060b8b9a99292505ef5d4becc3669dfa9d))

- Add argument to clean photometry that allows users to drop not host subtracted data
  ([`ea1fe31`](https://github.com/astro-otter/otter/commit/ea1fe314f00284516f05d42f7f5b7ea6b0aaea38))

- Allow users to pass additional kwargs to Transient.clean_photometry from Otter.get_phot
  ([`3ef876c`](https://github.com/astro-otter/otter/commit/3ef876cd8caab95716fda37ada0e647b24c0e0c6))

- Bug in filter name cleanup
  ([`6066926`](https://github.com/astro-otter/otter/commit/606692606090c3577d9046ba616caea697552580))

- Ensure the date_min/date_max cols are in the photometry df before filtering on them
  ([`dfc7229`](https://github.com/astro-otter/otter/commit/dfc7229b7ff7fea2756b55b4c4108f3f41b2690d))

- This is a patch for when date_min and date_max are null rather than NaN
  ([`83ec951`](https://github.com/astro-otter/otter/commit/83ec95134be32049158bc228b2e87366543d523d))


## v0.7.0 (2025-12-16)

### Bug Fixes

- Allow users to pass numpy arrays of names to Otter.query/Otter.get_meta (Issue #19)
  ([`ce89ccc`](https://github.com/astro-otter/otter/commit/ce89ccc3e13e6a9eb720199720b1998dfbd12b99))

- Fully address issue #32
  ([`63b7d0d`](https://github.com/astro-otter/otter/commit/63b7d0d2cf359b08040d59a879f021c2f2ba6421))

- Minor bug in photometry de-reddening
  ([`c6e5e33`](https://github.com/astro-otter/otter/commit/c6e5e33aa4718d158d32b9adf7a7320319bec0ec))

- Modernize test_transient unit tests
  ([`81e87e2`](https://github.com/astro-otter/otter/commit/81e87e25b1440f9d2e76f73516c71a6cad81b503))

- Plotter to use the prod db rather than dev
  ([`03351ce`](https://github.com/astro-otter/otter/commit/03351cec09ee8d09e71968d3fa30ef681a28052d))

- Redshift is now cast as a float (Issue #20)
  ([`7fba7a6`](https://github.com/astro-otter/otter/commit/7fba7a695e43d4fcf062cd914e385f94a0211e62))

### Features

- Add automatic MW dust extinction correction to Transient.clean_photometry
  ([`341fd2c`](https://github.com/astro-otter/otter/commit/341fd2c28e4b00d70b280d564a03f1aa7d2954d6))

- Standardize UVOIR filter names (Issue #21)
  ([`d798bbd`](https://github.com/astro-otter/otter/commit/d798bbde316316018cc6dbeade86f0dfe349cf0d))


## v0.6.2 (2025-12-12)

### Bug Fixes

- Add specific paths to only run CI unit tests when those files are edited
  ([`3b263ab`](https://github.com/astro-otter/otter/commit/3b263abe9167998614a84d8347a14b09c0aec32e))


## v0.6.1 (2025-12-12)

### Bug Fixes

- Ensure we don't run the semver workflow in an infinite loop
  ([`ec32d65`](https://github.com/astro-otter/otter/commit/ec32d6563ad5def4204892a92f5e3ac1a53d0f6b))


## v0.6.0 (2025-12-12)

### Bug Fixes

- A typo in the semver workflow
  ([`4666eaa`](https://github.com/astro-otter/otter/commit/4666eaa86fa4c9524daa0954c9f5df58301e0ad7))

- Add new semver.yml
  ([`cafbe26`](https://github.com/astro-otter/otter/commit/cafbe2648181ab067b4ed2976fd3931370d36033))

- Allow for versions to start with 0 with semantic-versioning
  ([`00bb113`](https://github.com/astro-otter/otter/commit/00bb113133bf358202b83e3e9f6ef04ee8bf5418))

- Following python-semantic-version docs add the token argument to actions/checkout
  ([`3a02c1f`](https://github.com/astro-otter/otter/commit/3a02c1f36e587f5134f827f4788a9e0b09815c68))

- Try adding in the repository name explicitly
  ([`43deea5`](https://github.com/astro-otter/otter/commit/43deea50cab952d2387c21f6b013a52059219c20))

- Try again with the semver workflow
  ([`6b7f948`](https://github.com/astro-otter/otter/commit/6b7f94866eb58a80f2fd17a6db5d3de6893beeb8))

- Try again with the semver workflow
  ([`c51e999`](https://github.com/astro-otter/otter/commit/c51e99971d768524b5d7786e05c1ecfec2c2c9f8))

- Try to fix the semver workflow file
  ([`9bb8bce`](https://github.com/astro-otter/otter/commit/9bb8bce1c642529f4c4658a14f1281d466b24dee))

- Try using a PAT instead of GITHUB_TOKEN
  ([`85deeff`](https://github.com/astro-otter/otter/commit/85deeffcd5fe638c865eab45cdc3f034cc1d177e))

- Update to use more recent versions of github actions
  ([`84926d8`](https://github.com/astro-otter/otter/commit/84926d87c5822b9a1163fb4fcb4fb58e06067b94))

### Features

- Try to add semantic versioning automatically
  ([`33f7664`](https://github.com/astro-otter/otter/commit/33f7664930b3742753fb67ad1906085f38515653))
