# Optional local site configuration

From the repository root:

```bash
cp examples/site/fastalake.site.example.sh fastalake.site.sh
# Edit fastalake.site.sh for your installation, then:
source fastalake.site.sh
```

The real `fastalake.site.sh` is ignored by Git. Keep local paths and credentials
out of the shared template. The manifest runner takes its data paths explicitly;
this template is optional. See [environment configuration](../../docs/get-started/environment.md).
