# The CeNGEN reference-data pins board

Returns a read-only `pins` board serving the `cengen-reference-data`
GitHub repo's published datasets over plain HTTPS (via
`raw.githubusercontent.com`), which hosts versioned CeNGEN reference
expression datasets (one per developmental stage / sex combination). New
datasets published to that repo become available to
[`list_cengen_datasets()`](https://cengenproject.github.io/cengenAnnotate/reference/list_cengen_datasets.md)
and
[`load_cengen_reference()`](https://cengenproject.github.io/cengenAnnotate/reference/load_cengen_reference.md)
with no package changes required.

## Usage

``` r
cengen_board(
  repo = getOption("cengenAnnotate.board_repo", "cengenproject/cengen-reference-data"),
  branch = getOption("cengenAnnotate.board_branch", "main"),
  token = Sys.getenv("GITHUB_PAT", unset = NA)
)
```

## Arguments

- repo:

  GitHub repo backing the board, as `"owner/repo"`.

- branch:

  Branch to read from.

- token:

  Optional GitHub personal access token, only needed if the repo is
  private. Defaults to the `GITHUB_PAT` environment variable;
  unauthenticated access is used if unset (fine for a public repo).

## Value

A pins board object.

## Details

This reads the repo's `_pins.yaml` manifest (see
[`pins::board_url()`](https://pins.rstudio.com/reference/board_url.html)),
which gives full dataset listing, metadata and versioning support
without needing git or a GitHub token for a public repo. Downloaded data
is cached locally by `pins` and only re-fetched when it changes.

Publishing new datasets is a separate, maintainer-side step (not part of
this package): clone `cengen-reference-data`, write pins to a local
[`pins::board_folder()`](https://pins.rstudio.com/reference/board_folder.html),
regenerate the manifest with
[`pins::write_board_manifest_yaml()`](https://pins.rstudio.com/reference/write_board_manifest_yaml.html),
and push — see that repo's README.
