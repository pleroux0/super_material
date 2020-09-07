#!/bin/sh

poetry run black --target-version $1 py38 super_material tests examples profile docs/conf.py
