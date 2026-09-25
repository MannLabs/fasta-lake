//! Capture provenance without mistaking a parent checkout for this source tree.
use std::path::Path;
use std::process::Command;
fn main() {
    println!("cargo:rerun-if-env-changed=FASTALAKE_BUILD_REV");
    println!("cargo:rerun-if-changed=src");
    let root = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap();
    let revision = std::env::var("FASTALAKE_BUILD_REV").unwrap_or_else(|_| {
        if !root.join(".git").exists() {
            return "source-archive".into();
        }
        let output = Command::new("git")
            .current_dir(root)
            .args(["rev-parse", "--short=12", "HEAD"])
            .output();
        let hash = output
            .ok()
            .filter(|o| o.status.success())
            .map(|o| String::from_utf8_lossy(&o.stdout).trim().to_string())
            .unwrap_or_else(|| "unknown".into());
        let dirty = Command::new("git")
            .current_dir(root)
            .args(["diff", "--quiet"])
            .status()
            .map(|s| !s.success())
            .unwrap_or(true);
        format!("{}{}", hash, if dirty { "-dirty" } else { "" })
    });
    println!("cargo:rustc-env=GIT_REV={}", revision);
}
