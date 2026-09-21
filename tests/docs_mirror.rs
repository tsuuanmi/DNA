use std::collections::BTreeSet;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

fn collect_files(root: &Path, extension: &str) -> io::Result<Vec<PathBuf>> {
    let mut files = Vec::new();
    let mut pending = vec![root.to_path_buf()];

    while let Some(directory) = pending.pop() {
        for entry in fs::read_dir(directory)? {
            let path = entry?.path();
            if path.is_dir() {
                pending.push(path);
            } else if path.extension().is_some_and(|value| value == extension) {
                files.push(path);
            }
        }
    }

    Ok(files)
}

fn relative_with_extension(path: &Path, root: &Path, extension: &str) -> io::Result<PathBuf> {
    let relative = path.strip_prefix(root).map_err(io::Error::other)?;
    let mut mapped = relative.to_path_buf();
    mapped.set_extension(extension);
    Ok(mapped)
}

#[test]
fn source_and_manual_docs_match_one_to_one() -> Result<(), Box<dyn std::error::Error>> {
    let repository = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let source_root = repository.join("src");
    let docs_root = repository.join("docs/src");

    let source_docs: BTreeSet<_> = collect_files(&source_root, "rs")?
        .into_iter()
        .map(|path| relative_with_extension(&path, &source_root, "md"))
        .collect::<Result<_, _>>()?;
    let manual_docs: BTreeSet<_> = collect_files(&docs_root, "md")?
        .into_iter()
        .map(|path| relative_with_extension(&path, &docs_root, "md"))
        .collect::<Result<_, _>>()?;

    assert_eq!(source_docs, manual_docs);
    Ok(())
}
