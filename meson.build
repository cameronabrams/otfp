project('otfp', 'c',
     meson_version : '>=0.57.0',
     version: run_command(meson.current_source_dir()/'version.sh', check: true).stdout().strip(),
)

meson.add_dist_script('version.sh')

subdir('src')

