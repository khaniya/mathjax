package com.example.mathjaxapplication;



import android.support.v7.app.ActionBarActivity;
import android.os.Bundle;
import android.view.Menu;
import android.view.MenuItem;
import android.os.*;
import android.view.*;
import android.widget.*;
import android.graphics.*;
import android.content.res.*;
import android.webkit.*;
import android.text.method.*;
import android.text.*;


public class MainActivity extends ActionBarActivity 
 implements View.OnClickListener{
	
	private String doubleEscapeTeX(String s) {
		String t="";
		for (int i=0; i < s.length(); i++) {
			if (s.charAt(i) == '\'') t += '\\';
			if (s.charAt(i) != '\n') t += s.charAt(i);
			if (s.charAt(i) == '\\') t += "\\";
		}
		return t;
	}
	
	private int exampleIndex = 0;
	private String getExample(int index) {
		return getResources().getStringArray(R.array.tex_examples)[index];
	}
	public void onClick(View v) {
		if (v == findViewById(R.id.button2)) {
			WebView w = (WebView) findViewById(R.id.webview);
			EditText e = (EditText) findViewById(R.id.edit);
			w.loadUrl("javascript:document.getElementById('math').innerHTML='\\\\["
			           +doubleEscapeTeX(e.getText().toString())+"\\\\]';");
			w.loadUrl("javascript:MathJax.Hub.Queue(['Typeset',MathJax.Hub]);");
		}
		else if (v == findViewById(R.id.button3)) {
			WebView w = (WebView) findViewById(R.id.webview);
			EditText e = (EditText) findViewById(R.id.edit);
			e.setText("");
			w.loadUrl("javascript:document.getElementById('math').innerHTML='';");
			w.loadUrl("javascript:MathJax.Hub.Queue(['Typeset',MathJax.Hub]);");
		}
		else if (v == findViewById(R.id.button4)) {
			WebView w = (WebView) findViewById(R.id.webview);
			EditText e = (EditText) findViewById(R.id.edit);
			e.setText(getExample(exampleIndex++));
			if (exampleIndex > getResources().getStringArray(R.array.tex_examples).length-1) 
				exampleIndex=0;
			w.loadUrl("javascript:document.getElementById('math').innerHTML='\\\\["
			          +doubleEscapeTeX(e.getText().toString())
					  +"\\\\]';");
			w.loadUrl("javascript:MathJax.Hub.Queue(['Typeset',MathJax.Hub]);");
		}
	}

    
	@Override
	protected void onCreate(Bundle savedInstanceState) {
		super.onCreate(savedInstanceState);
		setContentView(R.layout.activity_main);
		String[] a = getResources().getStringArray(R.array.tex_examples);
		WebView w = (WebView) findViewById(R.id.webview);
		w.getSettings().setJavaScriptEnabled(true);
		w.getSettings().setBuiltInZoomControls(true);
		w.loadDataWithBaseURL("http://bar", "<script type='text/x-mathjax-config'>"
                +"MathJax.Hub.Config({ " 
				  	+"showMathMenu: false, "
				  	+"jax: ['input/TeX','output/HTML-CSS'], "
				  	+"extensions: ['tex2jax.js'], " 
				  	+"TeX: { extensions: ['AMSmath.js','AMSsymbols.js',"
				  	  +"'noErrors.js','noUndefined.js'] } "
				  +"});</script>"
                +"<script type='text/javascript' "
				  +"src='file:///android_asset/MathJax/MathJax.js'"
				  +"></script><span id='math'></span>","text/html","utf-8","");
EditText e = (EditText) findViewById(R.id.edit);
e.setBackgroundColor(Color.LTGRAY);
e.setTextColor(Color.BLACK);
e.setText("");
Button b = (Button) findViewById(R.id.button2);
b.setOnClickListener(this);
b = (Button) findViewById(R.id.button3);
b.setOnClickListener(this);
b = (Button) findViewById(R.id.button4);
b.setOnClickListener(this);
//TextView t = (TextView) findViewById(R.id.textview3);
//t.setMovementMethod(LinkMovementMethod.getInstance());
//rt.setText(Html.fromHtml(t.getText().toString()));	

}

	

	@Override
	public boolean onCreateOptionsMenu(Menu menu) {
		// Inflate the menu; this adds items to the action bar if it is present.
		getMenuInflater().inflate(R.menu.main, menu);
		return true;
	}

	@Override
	public boolean onOptionsItemSelected(MenuItem item) {
		// Handle action bar item clicks here. The action bar will
		// automatically handle clicks on the Home/Up button, so long
		// as you specify a parent activity in AndroidManifest.xml.
		int id = item.getItemId();
		if (id == R.id.action_settings) {
			return true;
		}
		return super.onOptionsItemSelected(item);
	}
}
