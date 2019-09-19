
## List of version updates

### version-0.1.3
- Fix plotting issue in `compare_models`.
- `compare_models` can now directly take MLZ_model objects instead of a list, i.e., `compare_models(x, y)` instead of
`compare_models(list(x,y))`.
- Argument `type` in plot method for `MLZ_data` is added to produce one plot at a time.

### version-0.1.2
- Add an upper limit to Z as a function argument `Z.max`. Warning flag if model hits the upper limit.
- Brief description of MLmulti in documentation.

### version-0.1.1 
- Wrap examples that take longer than 5 secs to run with \dontrun{}. This is intended to address memory issues during CRAN check.
 
### version-0.1.0
- Initial release on CRAN.