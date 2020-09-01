#!/bin/sh

poetry run black --target-version py38 super_material tests examples profile
