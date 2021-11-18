from apero_tests.subtest import SubTest

class CustomBadpixTest(SubTest):
    def __init__(self, input_msg: str):
        """
        Subtest that does nothing, example for how we can define custom tests
        """
        self.msg = input_msg

        super().__init__(description=f"Test to show something")

    def run(self):

        self.result = self.msg + " :)"

        self.color = "Lime"

