# <b>EasySQ Overview</b>

### The purpose of EasySQ is to take a high level approach to Squidpy use.
<ul>
    <li>Get started with Squidpy in an easy to understand environment that can be easily translated into a more complex workflow.</li>
    <li>Abstract certain processes that are currently time consuming and complex in default Squidpy.</li>
    <li>Streamline and accelerate the process of outputing basic analysis.</li>
    <li>Work with more complex Squidpy processes through the use of both EasySQ and Squidpy.</li>
</ul>

# <b>EasySQ Installation and Setup</b>

This will walk you through how to install EasySQ from the github onto your IDE.<br>
If you are using PyCharm, follow the <i>PyCharm install guide.docx</i> located in the github repo. Otherwise, setup your IDE how you would like and then come back here.<br>
Once you have setup your IDE and run a "Hello World" script to make sure things are working, you are now ready to begin EasySQ installation.<br>

## Step 1 Getting from the GitHub.

We are going to download the repo as a .zip and place the files into the project manually.<br>
First, download the repo and unzip it.<br>
Secondly, place the following files into your project:<br>
<ol>
    <li><code>EasySQ.py</code></li>
    <li><code>EasySQTools.py</code></li>
    <li>Place the directory called <code>colors</code> into your project. This contains the color palettes for use in the demos. If you don't, you will need to make a directory and avoid using the provided palettes.</li>
    <li><code>requirements.txt</code> PyCharm will recognize this and ask to import the required packages. If your IDE does not, or PyCharm does not, please install these manually.</li>
</ol>

If this is your first time, I would also recommend:
<ol>
    <li><code>demo_1.py</code></li> for learning and working with EasySQ as it follows the Squidpy demos.
    <li><code>demo_2.py</code></li>
    <li><code>DisplayPalette.py</code></li> for showing the available color palettes.
</ol>

## Step 2 Running the Demos. (optional)

If you are new to EasySQ or Squidpy I would recommend looking through demos provided in the repo and by Squidpy. This will help things make sense when you try to make your own project.<br>
Start with <code>EasySQ Demo 1.ipynb</code> and then do <code>EasySQ Demo 2.ipynb</code>. These follow the linked Squidpy demos.
Once that is complete, I would recommend looking at the Squidpy Interface demo (once its completed), to see how to use EasySQ and Squidpy together.

## Step 3 Creating you own project.

Once you feel confident enough to create your own project, all you have to do is:<br>
<ol>
    <li><code>import EasySQ as esq</code></li>
    <li><code>if __name__ == "__main__":</code></li> EasySQ must be called from this if for some Squidpy functions to work.
    <li><code>esqAn = esq.Analysis(data_path='path/to/my/data/')</code></li> create a new Analysis class object and provide it a path to the directory for your data.
</ol>

And that's it! You can now start your EasySQ workflow.

### Demo 1 follows ["this"](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_vizgen.html) Squidpy tutorial on analyzing vizgen data.

I would recommend working through the provided "EasySQ Demo demo_number.ipynb." 1 and 2 are going to explain how to use EasySQ for your analysis of vizgen data.
There is also a demo on using EasySQ with Squidpy functions so that you can mix and match the two, which is also recommended.

### Demo 2 follows ["this"](https://squidpy.readthedocs.io/en/stable/notebooks/tutorials/tutorial_vizgen_mouse_liver.html) Squidpy tutorial on analyzing vizgen data.

#### For more information regarding any of the functions used, see the corresponding [Squidpy documentation](https://squidpy.readthedocs.io/en/stable/api.html#).

## <b>Disclaimer</b>

EasySQ is currently in an early development state.<br>
The tutorials and demos are the best way to learn EasySQ. Bugs may be present.<br>
For any questions, comments, or bugs, please contact Samuel Crouse at scrouse2@uwyo.edu.
