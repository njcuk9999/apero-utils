"""
DRS Tests that (try to) follow the APERO framework.

@author: vandalt
"""
import os
from datetime import datetime
from typing import Optional

from jinja2 import Environment, FileSystemLoader, select_autoescape

from . import OUTDIR, TEMPLATEDIR


def removext(name: str, ext: str = ".py"):
    """Remove extension of a file/recipe name

    :param name: file/recipe name
    :type name: str

    :returns: cleaned up name
    """
    if not ext.startswith("."):
        ext = "." + ext

    while name.endswith(ext):
        name = name[:-3]

    return name


class DrsTest:
    def __init__(
        self,
        instrument: Optional[str] = None,
        drs_recipe: Optional[str] = None,
        setup: Optional[str] = None,
    ):

        # get instrument
        self.instrument = instrument

        # Get setup path
        if setup is None:
            self.setup = os.environ["DRS_UCONFIG"]
        else:
            self.setup = setup

        # Set date at start of test
        self.date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")  # test date

        # Extract info from recipe
        if drs_recipe is not None:
            # Get recipe corresponding to test
            self.recipe = drs_recipe

            self.params = self.recipe.drs_params
            self.recipe_name = removext(self.recipe.name, ext=".py")
            self.ismaster = self.recipe.master

        else:
            self.recipe = "UnknownRecipe"
            self.params = None

        # TODO: Load log, outputs and calibdb in a nice way with APERO

    def gen_html(self, html_dict: dict):
        """Generate HTML summary from jinja2 template.

        :param html_dict: dict with all values used in html template
        :type html_dict: dict
        """

        # Jinja2 env
        env = Environment(
            loader=FileSystemLoader(TEMPLATEDIR),
            autoescape=select_autoescape(["html", "xml"]),
        )

        # Create template for test
        template = env.get_template(".".join([self.test_id, "html"]))

        html_text = template.render(html_dict)

        output_path = os.path.join(
            OUTDIR, self.test_id, ".".join([self.test_id, "html"])
        )
        with open(output_path, "w") as f:
            f.write(html_text)
