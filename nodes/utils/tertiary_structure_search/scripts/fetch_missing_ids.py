from pymol import cmd


def fetch_and_save(id, target_path):
    cmd.fetch(id)
    cmd.save(target_path, id)
